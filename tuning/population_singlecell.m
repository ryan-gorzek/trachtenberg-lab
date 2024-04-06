
%%%% add tools to path
addpath(genpath(fullfile("E:/Imaging/code/tools/")));

expst = [];
subj = "pv_vip_04";
sessions = ["110", "120", "130", "140", "150", "160"];
runs = ["000", "000", "000", "000", "000", "000"];
stims = ["trinoise", "randorisf", "battery1", "battery2", "battery3", "battery4"];
sid = 1:6;
stim_radius = 20; pix_per_deg = 14.3258;
eye_radius = 4;

for s = sid

    % load calcium data
    [rtdata, stdata, ops, stat] = load_calcium_data(subj, sessions(s), runs(s));
    
    % load stimulus data
    [frame_on, frame_off, params, stim_center] = load_stimulus_data(subj, sessions(s), runs(s));
    
    % load quadrature data
    quad = load_quad_data(subj, sessions(s), runs(s));
    
    % load eye data
    [eye, imgs] = load_eye_data(subj, sessions(s), runs(s));
    if stims(s) == "trinoise", eye_center = NaN; else, eye_center = expst.trinoise.spikes.eye.center; end
    
    % create experiment structure
    [expst.(stims(s)).spikes, expst.(stims(s)).dFF] = expstruct(stims(s), ...
                                                                rtdata, stdata, ops, stat, ...
                                                                frame_on, frame_off, params, ...
                                                                stim_center, pix_per_deg, ...
                                                                eye, imgs, eye_center, eye_radius, ...
                                                                quad);

    fprintf("Loaded %s...\n", stims(s));

end

%% make population summary

% ON domain criteria
on_dist = vertcat(expst.trinoise.spikes.stats.on_dist_fit) < 10;
on_snr = vertcat(expst.trinoise.spikes.stats.snr_on) > 5;
on_ve = vertcat(expst.trinoise.spikes.stats.cc_on) > 0.5;
on_centered = on_dist & on_snr & on_ve;
on_noise = ~(on_snr & on_ve);
% OFF domain criteria
off_dist = vertcat(expst.trinoise.spikes.stats.off_dist_fit) < 10;
off_snr = vertcat(expst.trinoise.spikes.stats.snr_off) > 5;
off_ve = vertcat(expst.trinoise.spikes.stats.cc_off) > 0.5;
off_centered = off_dist & off_snr & off_ve;
off_noise = ~(off_snr & off_ve);
% visually responsive
SNR = vertcat(expst.randorisf.spikes.stats.SNR) > expst.randorisf.spikes.calc.SNR_thr;
% intersect
cells = SNR & ((on_centered & off_centered) | (on_centered & off_noise) | (off_centered & on_noise));

f = figure; tiledlayout(6, 16);
% eye locations
nexttile(1, [2, 4]);
imshow(expst.trinoise.spikes.eye.img); hold on;
colors = hsv(6);
positions = []; labels = []; logicals = []; fractions = [];
for s = stims
    stimpos = expst.(s).spikes.eye.pos;
    positions = vertcat(positions, stimpos);
    labels = vertcat(labels, repmat(find(stims == s), [size(stimpos, 1), 1]));
    logicals = vertcat(logicals, expst.(s).spikes.eye.logical);
    fractions = vertcat(fractions, nnz(expst.(s).spikes.eye.logical) / numel(expst.(s).spikes.eye.logical));
end
positions_pass = positions(logicals == 1, :); labels_pass = labels(logicals == 1);
idx = randperm(size(positions_pass, 1));
scatter(positions_pass(idx, 1), positions_pass(idx, 2), 1, [0.9, 0.2, 0.2], "filled", "MarkerFaceAlpha", 0.2, "MarkerEdgeAlpha", 0.2);
scatter(positions(logicals == 0, 1), positions(logicals == 0, 2), 1, [0.2, 0.2, 0.9], "filled", "MarkerFaceAlpha", 0.2, "MarkerEdgeAlpha", 0.2);
% plot mask on top
start = expst.trinoise.spikes.eye.center - [eye_radius*2, eye_radius*2] ./ 2;
rectangle("Position", [start eye_radius*2 eye_radius*2], "Curvature", 0, "LineWidth", 1);
% plot stimuli modes
for s = stims
    stimmode = mode(expst.(s).spikes.eye.pos);
    textsize = ceil(fractions(stims == s) * 10);
    text(stimmode(1), stimmode(2), num2str(find(stims == s)), "FontSize", textsize, "HorizontalAlignment", "center", "VerticalAlignment", "middle");
end
% zoom in a bit
im_center = median(positions, 1);
range_x = 20; range_y = round((range_x / im_center(1)) * im_center(2));
xlim([im_center(1) - range_x, im_center(1) + range_x]); 
ylim([im_center(2) - range_y, im_center(2) + range_y]);

% imaging field
mask = 255 .* (sum(expst.trinoise.spikes.masks(:, :, cells), 3) > 0);
nmask = 255 .* (sum(expst.trinoise.spikes.masks(:, :, ~cells), 3) > 0);
mask = cat(3, mask, zeros(size(mask)), nmask);
nexttile(33, [2, 4]); 
imshow(repmat(expst.trinoise.spikes.maxproj(:, :, 1), [1, 1, 3])); hold on; h = imshow(mask); pbaspect([size(expst.trinoise.spikes.maxproj, [2, 1]), 1]); axis off; % 
set(h, "AlphaData", repmat(0.5, size(mask, [1, 2])));

% trinoise kernel locations
nexttile(65, [2, 4]);
imagesc(zeros(1080, 1920) + 0.2); colormap gray; hold on; pbaspect([1920, 1080, 1]); axis off;
ON = 80 .* [vertcat(expst.trinoise.spikes.stats.on_x_fit) + 1, vertcat(expst.trinoise.spikes.stats.on_y_fit) + 1];
OFF = 80 .* [vertcat(expst.trinoise.spikes.stats.off_x_fit) + 1, vertcat(expst.trinoise.spikes.stats.off_y_fit) + 1];
for c = 1:numel(cells)
    if cells(c) && on_centered(c)
        scatter(ON(c, 1), ON(c, 2), 10, "r", "filled", "MarkerFaceAlpha", 0.5, "MarkerEdgeAlpha", 0.5);
    elseif cells(c) && off_centered(c)
        scatter(OFF(c, 1), OFF(c, 2), 10, "r", "filled", "MarkerFaceAlpha", 0.2, "MarkerEdgeAlpha", 0.2);
    elseif ~cells(c)
        scatter(ON(c, 1), ON(c, 2), 10, "b", "filled", "MarkerFaceAlpha", 0.2, "MarkerEdgeAlpha", 0.2);
    end
end
% color by grid location
% grid_colors = hot(9);
% grid = vertcat(expst.randorisf.spikes.cells.grid);
% for c = 1:numel(cells)
%     if cells(c), scatter(ON(c, 1), ON(c, 2), 10, "MarkerFaceColor", grid_colors(grid(c), :), "MarkerEdgeColor", grid_colors(grid(c), :), ...
%                                                                   "MarkerFaceAlpha", 0.5, "MarkerEdgeAlpha", 0.5); end
%     if cells(c), scatter(OFF(c, 1), OFF(c, 2), 10, "MarkerFaceColor", grid_colors(grid(c), :), "MarkerEdgeColor", grid_colors(grid(c), :), ...
%                                                                       "MarkerFaceAlpha", 0.5, "MarkerEdgeAlpha", 0.5); end
% end
% plot stimulus location on top
pixel_radius = stim_radius * pix_per_deg;
start = expst.trinoise.spikes.stim_center - [pixel_radius, pixel_radius] ./ 2;
rectangle("Position", [start pixel_radius pixel_radius], "Curvature", 1);
rectangle("Position", [1, 1, [1920, 1080] - 1], "Curvature", 0);

% randorisf histograms incl. ori, sf, SNR (and thr), F1F0
nexttile(5, [2, 2]); hold on; feature_hist(expst.randorisf.spikes.table, "SNR_ROSF_spikes", true(size(cells)), "SNR"); xline(expst.randorisf.spikes.calc.SNR_thr, "--r", "LineWidth", 1);
nexttile(37, [2, 2]); hold on; feature_hist(expst.randorisf.spikes.table, "F1F0_ROSF_spikes", cells, "F1F0");
nexttile(69, [2, 2]); hold on; feature_hist(expst.randorisf.spikes.table, "Ori_ROSF_spikes", cells, "Orientation");
nexttile(7, [2, 2]); hold on; feature_hist(expst.randorisf.spikes.table, "SF_ROSF_spikes", cells, "Spatial Frequency");

% battery1 histograms, incl. dir, sf, DSI
nexttile(39, [2, 2]); hold on; feature_hist(expst.battery1.spikes.table, "Dir_B1_spikes", cells, "Direction");
nexttile(71, [2, 2]); hold on; feature_hist(expst.battery1.spikes.table, "SF_B1_spikes", cells, "Spatial Frequency");
nexttile(9, [2, 2]); hold on; feature_hist(expst.battery1.spikes.table, "DSI_B1_spikes", cells, "DSI");

% battery2 histograms, incl. tf
nexttile(41, [2, 2]); hold on; feature_hist(expst.battery2.spikes.table, "TF_B2_spikes", cells, "Temporal Frequency");

% battery3 histograms, incl. contrast, size
nexttile(73, [2, 2]); hold on; feature_hist(expst.battery3.spikes.table, "Ctrst_B3_spikes", cells, "Contrast");
nexttile(11, [2, 2]); hold on; feature_hist(expst.battery3.spikes.table, "Size_B3_spikes", cells, "Size");
% running/stationary size tuning curve
nexttile(43, [2, 2]); axis square; axis tight; hold on;
stationary = mean(vertcat(expst.battery3.dFF.stats(cells).size_stationary), 1, "omitnan"); running = mean(vertcat(expst.battery3.dFF.stats(cells).size_running), 1, "omitnan");
plot(expst.battery3.spikes.viz.ticks.cont{3}, stationary, "Color", [0.7, 0.7, 0.7]);
scatter(expst.battery3.spikes.viz.ticks.cont{3}, stationary, "filled", "MarkerFaceColor", [1, 1, 1], "MarkerEdgeColor", [0.7, 0.7, 0.7]);
plot(expst.battery3.spikes.viz.ticks.cont{3}, running, "Color", [0.7, 0.7, 0.7]);
scatter(expst.battery3.spikes.viz.ticks.cont{3}, running, "filled", "MarkerFaceColor", [0.7, 0.7, 0.7], "MarkerEdgeColor", [0.7, 0.7, 0.7]);
xticks(expst.battery3.spikes.viz.ticks.cont{3}); xticklabels(expst.battery3.spikes.viz.ticklabels{3}); xlabel(expst.battery3.spikes.viz.labels{3});
ylabel("\DeltaF/F");
legend("", "Stationary", "", "Running", "Box", "off", "Location", "northoutside", "Orientation", "horizontal");

% cross vs iso
X = vertcat(expst.battery4.dFF.stats(cells).cross_mean);
Y = vertcat(expst.battery4.dFF.stats(cells).iso_all_mean);
nexttile(13, [2, 2]); hold on;
xline(0, "--k", "LineWidth", 1);
yline(0, "--k", "LineWidth", 1);
plot(-1:1, -1:1, "--k", "LineWidth", 1);
X(X > 0.5) = 0.5;
Y(Y > 0.5) = 0.5;
scatter(X(:, 3), Y(:, 3), 10, [0, 0, 0], "filled", "MarkerFaceAlpha", 0.3, "MarkerEdgeAlpha", 0.3);
xlim([-0.1, 0.5]); xticks(-0.1:0.1:0.5);
ylim([-0.1, 0.5]); yticks(-0.1:0.1:0.5);
xlabel("Cross [\DeltaF/F]"); ylabel("Iso [\DeltaF/F]");
axis square;

% center vs iso
X = vertcat(expst.battery4.dFF.stats(cells).center_mean);
Y = vertcat(expst.battery4.dFF.stats(cells).iso_all_mean);
nexttile(15, [2, 2]); hold on;
xline(0, "--k", "LineWidth", 1);
yline(0, "--k", "LineWidth", 1);
plot(-1:1, -1:1, "--k", "LineWidth", 1);
X(X > 0.5) = 0.5;
Y(Y > 0.5) = 0.5;
scatter(X(:, 3), Y(:, 3), 10, [0, 0, 0], "filled", "MarkerFaceAlpha", 0.3, "MarkerEdgeAlpha", 0.3);
xlim([-0.1, 0.5]); xticks(-0.1:0.1:0.5);
ylim([-0.1, 0.5]); yticks(-0.1:0.1:0.5);
xlabel("Center [\DeltaF/F]"); ylabel("Iso [\DeltaF/F]");
axis square;

% cross/iso tuning curves at different sizes
nexttile(77, [2, 2]); axis square; axis tight; hold on;
iso = mean(vertcat(expst.battery4.dFF.stats(cells).iso_mean), 1); cross = mean(vertcat(expst.battery4.dFF.stats(cells).cross_mean), 1);
plot(expst.battery4.spikes.viz.ticks.cont{5}, iso, "Color", [0.7, 0.7, 0.7]);
scatter(expst.battery4.spikes.viz.ticks.cont{5}, iso, "filled", "MarkerFaceColor", [1, 1, 1], "MarkerEdgeColor", [0.7, 0.7, 0.7]);
plot(expst.battery4.spikes.viz.ticks.cont{5}, cross, "Color", [0.7, 0.7, 0.7]);
scatter(expst.battery4.spikes.viz.ticks.cont{5}, cross, "filled", "MarkerFaceColor", [0.7, 0.7, 0.7], "MarkerEdgeColor", [0.7, 0.7, 0.7]);
xticks(expst.battery4.dFF.viz.ticks.cont{5}); xticklabels(expst.battery4.dFF.viz.ticklabels{5}); xlabel(expst.battery4.dFF.viz.labels{5});
ylabel("\DeltaF/F");
legend("", "Iso", "", "Cross", "Box", "off", "Location", "northoutside", "Orientation", "horizontal");
% center/surround tuning curves at different sizes
nexttile(75, [2, 2]); axis square; axis tight; hold on;
surround = mean(vertcat(expst.battery4.dFF.stats(cells).surround_mean), 1); center = mean(vertcat(expst.battery4.dFF.stats(cells).center_mean), 1);
plot(expst.battery4.spikes.viz.ticks.cont{5}, surround, "Color", [0.7, 0.7, 0.7]);
scatter(expst.battery4.spikes.viz.ticks.cont{5}, surround, "filled", "MarkerFaceColor", [1, 1, 1], "MarkerEdgeColor", [0.7, 0.7, 0.7]);
plot(expst.battery4.spikes.viz.ticks.cont{5}, center, "Color", [0.7, 0.7, 0.7]);
scatter(expst.battery4.spikes.viz.ticks.cont{5}, center, "filled", "MarkerFaceColor", [0.7, 0.7, 0.7], "MarkerEdgeColor", [0.7, 0.7, 0.7]);
xticks(expst.battery4.dFF.viz.ticks.cont{5}); xticklabels(expst.battery4.dFF.viz.ticklabels{5}); xlabel(expst.battery4.dFF.viz.labels{5});
ylabel("\DeltaF/F");
legend("", "Surround", "", "Center", "Box", "off", "Location", "northoutside", "Orientation", "horizontal");

% battery4 histograms, incl. modulation index, ITI
nexttile(45, [2, 2]); hold on; feature_hist(expst.battery4.dFF.table, "cri_mod_mean_B4_dFF", cells, "Cross/Iso Modulation Index");
nexttile(47, [2, 2]); hold on; feature_hist(expst.battery4.dFF.table, "ci_mod_mean_B4_dFF", cells, "Center/Iso Modulation Index");

set(f, "Color", "w");
f.Position = [219 347 1653 704];

%% make cell summary

[~, idx] = sort(vertcat(expst.trinoise.spikes.stats.on_dist_fit) + vertcat(expst.trinoise.spikes.stats.off_dist_fit), "ascend");
idx = idx(ismember(idx, find(cells)));
id = idx(1);

f = figure; tiledlayout(4, 15);
% trinoise kernel
nexttile(1, [2, 3]);
RdBu_r = colorMap({[0, 0, 1], [1, 1, 1], [1, 0, 0]});
% ON
kern_on = repelem(expst.trinoise.spikes.stats(id).kern_on_fit, 80, 80); % 
xpad_on = ones(size(kern_on, 1), (1920 - size(kern_on, 2)) / 2) .* kern_on(1, 1);
kern_on = horzcat(xpad_on, kern_on, xpad_on);
ypad_top_on = ones((1080 - size(kern_on, 1)) / 2.5, size(kern_on, 2)) .* kern_on(1, 1);
ypad_bot_on = ones(1.5 * ((1080 - size(kern_on, 1)) / 2.5), size(kern_on, 2)) .* kern_on(1, 1);
kern_on = vertcat(ypad_top_on, kern_on, ypad_bot_on);
% OFF
kern_off = repelem(expst.trinoise.spikes.stats(id).kern_off_fit, 80, 80); % 
xpad_off = ones(size(kern_off, 1), (1920 - size(kern_off, 2)) / 2) .* kern_off(1, 1);
kern_off = horzcat(xpad_off, kern_off, xpad_off);
ypad_top_off = ones((1080 - size(kern_off, 1)) / 2.5, size(kern_off, 2)) .* kern_off(1, 1);
ypad_bot_off = ones(1.5 * ((1080 - size(kern_off, 1)) / 2.5), size(kern_off, 2)) .* kern_off(1, 1);
kern_off = vertcat(ypad_top_off, kern_off, ypad_bot_off);
% plot
limit = min([max(kern_off, [], "all"), max(kern_on, [], "all")]);
imagesc(kern_on - kern_off); hold on; colormap(RdBu_r); pbaspect([flip(size(kern_on)), 1]); axis off; colorbar("westoutside"); clim([-limit, limit]);
ON = 80 .* [expst.trinoise.spikes.stats(id).on_x_fit + 1, expst.trinoise.spikes.stats(id).on_y_fit + 1]; %
OFF = 80 .* [expst.trinoise.spikes.stats(id).off_x_fit + 1, expst.trinoise.spikes.stats(id).off_y_fit + 1]; %
yticklabels([]);
scatter(ON(1), ON(2), "+k"); scatter(OFF(1), OFF(2), "+k");
% plot stimulus location on top
pixel_radius = stim_radius * pix_per_deg;
start = expst.trinoise.spikes.stim_center - [pixel_radius, pixel_radius] ./ 2;
rectangle("Position", [start pixel_radius pixel_radius], "Curvature", 1);
rectangle("Position", [1, 1, flip(size(kern_on) - 1)], "Curvature", 0);

% randorisf kernel
nexttile(31, [2, 3]);
imagesc(expst.randorisf.spikes.kernels(id).kernsmooth'); pbaspect([size(expst.randorisf.spikes.kernels(id).kern), 1]); colormap(gca, "parula");
xticks(expst.randorisf.spikes.viz.ticks.disc{1}); xticklabels(expst.randorisf.spikes.viz.ticklabels{1}); xlabel(expst.randorisf.spikes.viz.labels{1});
yticks(expst.randorisf.spikes.viz.ticks.disc{2}); yticklabels(expst.randorisf.spikes.viz.ticklabels{2}); ylabel(expst.randorisf.spikes.viz.labels{2});

% battery1 kernel
nexttile(4, [2, 2]);
imagesc(expst.battery1.spikes.kernels(id).kernsmooth); pbaspect([flip(size(expst.battery1.spikes.kernels(id).kern)), 1]); colormap(gca, "parula");
xticks(expst.battery1.spikes.viz.ticks.disc{2}); xticklabels(expst.battery1.spikes.viz.ticklabels{2}); xlabel(expst.battery1.spikes.viz.labels{2});
yticks(expst.battery1.spikes.viz.ticks.disc{1}); yticklabels(expst.battery1.spikes.viz.ticklabels{1}); ylabel(expst.battery1.spikes.viz.labels{1});

% battery2 kernel
nexttile(34, [2, 2]);
imagesc(expst.battery2.spikes.kernels(id).kernsmooth); pbaspect([flip(size(expst.battery2.spikes.kernels(id).kern)), 1]); colormap(gca, "parula");
xticks(expst.battery2.spikes.viz.ticks.disc{2}); xticklabels(expst.battery2.spikes.viz.ticklabels{2}); xlabel(expst.battery2.spikes.viz.labels{2});
yticks(expst.battery2.spikes.viz.ticks.disc{1}); yticklabels(expst.battery2.spikes.viz.ticklabels{1}); ylabel(expst.battery2.spikes.viz.labels{1});

% battery3
% contrast/size kernel
nexttile(6, [2, 2]);
imagesc(squeeze(expst.battery3.spikes.kernels(id).kernsmooth(:, expst.battery3.spikes.kernels(id).peakidx(2), :))); colormap(gca, "parula");
pbaspect([flip(size(expst.battery3.spikes.kernels(id).kern, [1, 3])), 1]);
xticks(expst.battery3.spikes.viz.ticks.disc{3}); xticklabels(expst.battery3.spikes.viz.ticklabels{3}); xlabel(expst.battery3.spikes.viz.labels{3});
yticks(expst.battery3.spikes.viz.ticks.disc{1}); yticklabels(expst.battery3.spikes.viz.ticklabels{1}); ylabel(expst.battery3.spikes.viz.labels{1});
% contrast/ori kernel
nexttile(8, [2, 2]);
imagesc(squeeze(expst.battery3.spikes.kernels(id).kernsmooth(:, :, expst.battery3.spikes.kernels(id).peakidx(3)))); colormap(gca, "parula");
pbaspect([flip(size(expst.battery3.spikes.kernels(id).kern, [1, 2])), 1]);
xticks(expst.battery3.spikes.viz.ticks.disc{2}); xticklabels(expst.battery3.spikes.viz.ticklabels{2}); xlabel(expst.battery3.spikes.viz.labels{2});
yticks(expst.battery3.spikes.viz.ticks.disc{1}); yticklabels(expst.battery3.spikes.viz.ticklabels{1}); ylabel(expst.battery3.spikes.viz.labels{1});
% size/ori kernel
nexttile(36, [2, 2]);
imagesc(squeeze(expst.battery3.spikes.kernels(id).kernsmooth(expst.battery3.spikes.kernels(id).peakidx(1), :, :))'); colormap(gca, "parula");
pbaspect([flip(size(expst.battery3.spikes.kernels(id).kern, [3, 2])), 1]);
xticks(expst.battery3.spikes.viz.ticks.disc{2}); xticklabels(expst.battery3.spikes.viz.ticklabels{2}); xlabel(expst.battery3.spikes.viz.labels{2});
yticks(expst.battery3.spikes.viz.ticks.disc{3}); yticklabels(expst.battery3.spikes.viz.ticklabels{3}); ylabel(expst.battery3.spikes.viz.labels{3});
% running/stationary size tuning curve
nexttile(38, [2, 2]); axis square; axis tight; hold on;
stationary = expst.battery3.dFF.stats(id).size_stationary; running = expst.battery3.dFF.stats(id).size_running;
plot(expst.battery3.spikes.viz.ticks.cont{3}, stationary, "Color", [0.7, 0.7, 0.7]);
scatter(expst.battery3.spikes.viz.ticks.cont{3}, stationary, "filled", "MarkerFaceColor", [1, 1, 1], "MarkerEdgeColor", [0.7, 0.7, 0.7]);
plot(expst.battery3.spikes.viz.ticks.cont{3}, running, "Color", [0.7, 0.7, 0.7]);
scatter(expst.battery3.spikes.viz.ticks.cont{3}, running, "filled", "MarkerFaceColor", [0.7, 0.7, 0.7], "MarkerEdgeColor", [0.7, 0.7, 0.7]);
xticks(expst.battery3.spikes.viz.ticks.cont{3}); xticklabels(expst.battery3.spikes.viz.ticklabels{3}); xlabel(expst.battery3.spikes.viz.labels{3});
ylabel("\DeltaF/F");
legend("", "Stationary", "", "Running", "Box", "off", "Location", "northoutside", "Orientation", "horizontal");

% context modulation
sz = 20;
nexttile(10, [2, 2]); hold on;
ylims = ctxmod_psth(expst, id, "Center", sz, NaN);
nexttile(12, [2, 2]); hold on;
ctxmod_psth(expst, id, "Surround", sz, ylims);
nexttile(40, [2, 2]); hold on;
ctxmod_psth(expst, id, "Cross", sz, ylims);
nexttile(42, [2, 2]); hold on;
ctxmod_psth(expst, id, "Iso", sz, ylims);
% center/surround tuning curves at different sizes
nexttile(14, [2, 2]); axis square; axis tight; hold on;
surround = expst.battery4.dFF.stats(id).surround_mean; center = expst.battery4.dFF.stats(id).center_mean; 
plot(expst.battery4.spikes.viz.ticks.cont{5}, surround, "Color", [0.7, 0.7, 0.7]);
scatter(expst.battery4.spikes.viz.ticks.cont{5}, surround, "filled", "MarkerFaceColor", [1, 1, 1], "MarkerEdgeColor", [0.7, 0.7, 0.7]);
plot(expst.battery4.spikes.viz.ticks.cont{5}, center, "Color", [0.7, 0.7, 0.7]);
scatter(expst.battery4.spikes.viz.ticks.cont{5}, center, "filled", "MarkerFaceColor", [0.7, 0.7, 0.7], "MarkerEdgeColor", [0.7, 0.7, 0.7]);
xticks(expst.battery4.dFF.viz.ticks.cont{5}); xticklabels(expst.battery4.dFF.viz.ticklabels{5}); xlabel(expst.battery4.dFF.viz.labels{5});
ylabel("\DeltaF/F");
legend("", "Surround", "", "Center", "Box", "off", "Location", "northoutside", "Orientation", "horizontal");
% cross/iso tuning curves at different sizes
nexttile(44, [2, 2]); axis square; axis tight; hold on;
iso = expst.battery4.dFF.stats(id).iso_mean; cross = expst.battery4.dFF.stats(id).cross_mean;
plot(expst.battery4.spikes.viz.ticks.cont{5}, iso, "Color", [0.7, 0.7, 0.7]);
scatter(expst.battery4.spikes.viz.ticks.cont{5}, iso, "filled", "MarkerFaceColor", [1, 1, 1], "MarkerEdgeColor", [0.7, 0.7, 0.7]);
plot(expst.battery4.spikes.viz.ticks.cont{5}, cross, "Color", [0.7, 0.7, 0.7]);
scatter(expst.battery4.spikes.viz.ticks.cont{5}, cross, "filled", "MarkerFaceColor", [0.7, 0.7, 0.7], "MarkerEdgeColor", [0.7, 0.7, 0.7]);
xticks(expst.battery4.dFF.viz.ticks.cont{5}); xticklabels(expst.battery4.dFF.viz.ticklabels{5}); xlabel(expst.battery4.dFF.viz.labels{5});
ylabel("\DeltaF/F");
legend("", "Iso", "", "Cross", "Box", "off", "Location", "northoutside", "Orientation", "horizontal");

set(f, "Color", "w");
f.Position = [121 614.7143 1.4314e+03 458.8571];

%% check on stationary/running size tuning

id = idx(1);

PSTHs = expst.battery3.dFF.psths(id).mat;
stims = expst.battery3.dFF.stimdata.param(expst.battery3.dFF.stimdata.idx, :);
stimidx = expst.battery3.dFF.stimdata.stimidx;
runidx = vertcat(expst.battery3.dFF.stims.run);

stationary = PSTHs(runidx == 0, :);
somerun = PSTHs(runidx == 1, :);
running = PSTHs(runidx == 2, :);

figure; tiledlayout(1, 3);
nexttile; plot(stationary');
nexttile; plot(somerun');
nexttile; plot(running');

%% plot PSTHs

id = 1;

f = figure; tiledlayout(2, 4);
nexttile; hold on; psth(expst.battery1.dFF, id);
nexttile; hold on; psth(expst.battery2.dFF, id);
nexttile; hold on; psth(expst.battery3.dFF, id);
nexttile; hold on; psth(expst.battery4.dFF, id);
nexttile; hold on; psth(expst.battery1.spikes, id);
nexttile; hold on; psth(expst.battery2.spikes, id);
nexttile; hold on; psth(expst.battery3.spikes, id);
nexttile; hold on; psth(expst.battery4.spikes, id);

%% play movie of eye and tracking

figure;
for fr = 1:size(imgs, 3)
    imshow(imgs(:, :, fr)); hold on; scatter(eye(fr).Centroid(1), eye(fr).Centroid(2), "+r"); axis off;
    pause(0.001);
    drawnow;
end

%% functions

function psth(strct, id)

    data = strct.psths(id).mat;
    
    yline(0, "--r", "LineWidth", 1);
    xline(-15, "--k", "LineWidth", 1);
    xline(0, "--k", "LineWidth", 1);
    xline(15, "--k", "LineWidth", 1);
    wdw = strct.psths(id).window(1) : strct.psths(id).window(2);
    stdev = std(data, 0, 1);
    fill([wdw, flip(wdw)], [mean(data, 1) + stdev, flip(mean(data, 1)) - flip(stdev)], [0.7, 0.7, 0.7]);
    plot(strct.psths(id).window(1) : strct.psths(id).window(2), mean(data, 1), "LineWidth", 2, "Color", "k");
    axis square; axis tight;
    xlim([-5, strct.psths(id).window(2)]);
    ylim([-0.25, 0.5]);

end

function ylims = ctxmod_psth(expst, id, context, sz, ylims)

    sizes = expst.battery4.dFF.stimdata.param(expst.battery4.dFF.stimdata.idx, 5);
    data = expst.battery4.dFF.psths(id).mat(expst.battery4.dFF.stimdata.context.(lower(context)) & sizes == sz, :);
    
    yline(0, "--r", "LineWidth", 1);
    xline(-15, "--k", "LineWidth", 1);
    xline(0, "--k", "LineWidth", 1);
    xline(15, "--k", "LineWidth", 1);
    wdw = expst.battery4.dFF.psths(id).window(1) : expst.battery4.dFF.psths(id).window(2);
    stdev = std(data, 0, 1);
    fill([wdw, flip(wdw)], [mean(data, 1) + stdev, flip(mean(data, 1)) - flip(stdev)], [0.7, 0.7, 0.7]);
    plot(expst.battery4.dFF.psths(id).window(1) : expst.battery4.dFF.psths(id).window(2), mean(data, 1), "LineWidth", 2, "Color", "k");
    title(context);
    axis square; axis tight;
    xlim([-5, 30]);
    if ~isnan(ylims), ylim(ylims); end
    ylims = gca().YLim;
    ylabel("\DeltaF/F");

end

function feature_hist(tbl, featname, idx, t)
    feat = ismember(tbl.Properties.VariableNames, featname);
    scale = string(tbl.Properties.CustomProperties.FeatureScale{feat});
    binedges = tbl.Properties.CustomProperties.BinEdges{feat};
    ticks = tbl.Properties.CustomProperties.Ticks{feat};
    ticklabels = tbl.Properties.CustomProperties.TickLabels{feat};
    if scale == "linear" || scale == "polar"
        histogram(tbl.(featname)(idx), binedges, "Normalization", "probability");
    elseif scale == "log"
        histogram(log10(tbl.(featname)(idx)), binedges, "Normalization", "probability");
    end
    xticks(ticks); xticklabels(ticklabels);
    title(t, "Interpreter", "none");
    axis square; axis tight;
end
