
%%%% add tools to path
addpath(genpath(fullfile("E:/Imaging/code/tools/")));

expst = [];
subj = "pv_vip_01";
sessions = ["410", "420", "430", "440", "450", "460"];
runs = ["000", "000", "000", "000", "000", "000"];
stims = ["trinoise", "randorisf", "battery1", "battery2", "battery3", "battery4"];
sid = 1:6;
stim_center = [960, 540]; stim_radius = 22; pix_per_deg = 14.3258;

for s = sid

    % load calcium data
    [rtdata, stdata, xy, grid] = load_calcium_data(subj, sessions(s), runs(s));
    
    % load stimulus data
    [frame_on, frame_off, params] = load_stimulus_data(subj, sessions(s), runs(s));
    
    % load quadrature data
    quad = load_quad_data(subj, sessions(s), runs(s));
    
    % load eye data
    eye = load_eye_data(subj, sessions(s), runs(s));
    
    % create experiment structure
    [expst.(stims(s)).spikes, expst.(stims(s)).dFF] = expstruct(stims(s), rtdata, stdata, xy, grid, frame_on, frame_off, params, eye, quad);

end

%% make cell summary

[~, idx] = sort(vertcat(expst.randorisf.spikes.stats.SNR), "descend");
id = idx(60);

f = figure; tiledlayout(4, 13);
% trinoise kernel
nexttile(1, [2, 3]);
RdBu_r = colorMap({[0, 0, 1], [1, 1, 1], [1, 0, 0]});
kern_on = repelem(expst.trinoise.spikes.stats(id).kern_on_fit, 80, 80);
kern_on = horzcat(kern_on, ones(size(kern_on, 1), 1920 - size(kern_on, 2)) .* kern_on(1, 1));
kern_on = vertcat(kern_on, ones(1080 - size(kern_on, 1), size(kern_on, 2)) .* kern_on(1, 1));
kern_off = repelem(expst.trinoise.spikes.stats(id).kern_off_fit, 80, 80);
kern_off = horzcat(kern_off, ones(size(kern_off, 1), 1920 - size(kern_off, 2)) .* kern_off(1, 1));
kern_off = vertcat(kern_off, ones(1080 - size(kern_off, 1), size(kern_off, 2)) .* kern_off(1, 1));
limit = min([max(kern_off, [], "all"), max(kern_on, [], "all")]);
imagesc(kern_on - kern_off); hold on; colormap(RdBu_r); pbaspect([flip(size(kern_on)), 1]); axis off; colorbar("westoutside"); clim([-limit, limit]);
ON = 80 .* [expst.trinoise.spikes.stats(id).on_x_fit, expst.trinoise.spikes.stats(id).on_y_fit] - 40;
OFF = 80 .* [expst.trinoise.spikes.stats(id).off_x_fit, expst.trinoise.spikes.stats(id).off_y_fit] - 40;
yticklabels([]);
scatter(ON(1), ON(2), "+k"); scatter(OFF(1), OFF(2), "+k");
% plot stimulus location on top
pixel_radius = stim_radius * pix_per_deg;
start = stim_center - [pixel_radius, pixel_radius] ./ 2;
rectangle("Position", [start pixel_radius pixel_radius], "Curvature", 1);
rectangle("Position", [1, 1, flip(size(kern_on) - 1)], "Curvature", 0);

% randorisf kernel
nexttile(27, [2, 3]);
imagesc(expst.randorisf.spikes.kernels(id).kernsmooth'); pbaspect([size(expst.randorisf.spikes.kernels(id).kern), 1]); colormap(gca, "parula");
xticks(expst.randorisf.spikes.viz.ticks.disc{1}); xticklabels(expst.randorisf.spikes.viz.ticklabels{1}); xlabel(expst.randorisf.spikes.viz.labels{1});
yticks(expst.randorisf.spikes.viz.ticks.disc{2}); yticklabels(expst.randorisf.spikes.viz.ticklabels{2}); ylabel(expst.randorisf.spikes.viz.labels{2});

% battery1 kernel
nexttile(4, [2, 2]);
imagesc(expst.battery1.spikes.kernels(id).kernsmooth); pbaspect([flip(size(expst.battery1.spikes.kernels(id).kern)), 1]); colormap(gca, "parula");
xticks(expst.battery1.spikes.viz.ticks.disc{2}); xticklabels(expst.battery1.spikes.viz.ticklabels{2}); xlabel(expst.battery1.spikes.viz.labels{2});
yticks(expst.battery1.spikes.viz.ticks.disc{1}); yticklabels(expst.battery1.spikes.viz.ticklabels{1}); ylabel(expst.battery1.spikes.viz.labels{1});

% battery2 kernel
nexttile(30, [2, 2]);
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
nexttile(32, [2, 2]);
imagesc(squeeze(expst.battery3.spikes.kernels(id).kernsmooth(expst.battery3.spikes.kernels(id).peakidx(1), :, :))'); colormap(gca, "parula");
pbaspect([flip(size(expst.battery3.spikes.kernels(id).kern, [3, 2])), 1]);
xticks(expst.battery3.spikes.viz.ticks.disc{2}); xticklabels(expst.battery3.spikes.viz.ticklabels{2}); xlabel(expst.battery3.spikes.viz.labels{2});
yticks(expst.battery3.spikes.viz.ticks.disc{3}); yticklabels(expst.battery3.spikes.viz.ticklabels{3}); ylabel(expst.battery3.spikes.viz.labels{3});
% running/stationary size tuning curve
nexttile(34, [2, 2]); axis square; axis tight; hold on;
stationary = expst.battery3.spikes.stats(id).size_stationary; running = expst.battery3.spikes.stats(id).size_running;
plot(expst.battery3.spikes.viz.ticks.cont{3}, stationary, "Color", [0.7, 0.7, 0.7]);
scatter(expst.battery3.spikes.viz.ticks.cont{3}, stationary, "filled", "MarkerFaceColor", [1, 1, 1], "MarkerEdgeColor", [0.7, 0.7, 0.7]);
plot(expst.battery3.spikes.viz.ticks.cont{3}, running, "Color", [0.7, 0.7, 0.7]);
scatter(expst.battery3.spikes.viz.ticks.cont{3}, running, "filled", "MarkerFaceColor", [0.7, 0.7, 0.7], "MarkerEdgeColor", [0.7, 0.7, 0.7]);
xticks(expst.battery3.spikes.viz.ticks.cont{3}); xticklabels(expst.battery3.spikes.viz.ticklabels{3}); xlabel(expst.battery3.spikes.viz.labels{3});
ylabel("Spike Probability");
legend("", "Stationary", "", "Running", "Box", "off", "Location", "northoutside", "Orientation", "horizontal");

% context modulation
nexttile(9, [2, 2]); hold on;

set(f, "Color", "w");
f.Position = [121 601.5714 1072 472.0000];

function ctxmod_psth(id, context, sz)

    sizes = expst.battery4.dFF.stimdata.param(activity.stimdata.idx, 5);
    data = expst.battery4.dFF.psths(id).mat(expst.battery4.dFF.stimdata.context.(lower(context)) & sizes == sz, :);
    
    yline(0, "--r", "LineWidth", 1);
    xline(-15, "--k", "LineWidth", 2);
    xline(0, "--k", "LineWidth", 2);
    xline(15, "--k", "LineWidth", 2);
    plot(activity.psths(1).window(1) : activity.psths(1).window(2), data.(st), "LineWidth", 1, "Color", [0.8, 0.8, 0.8, 0.7]);
    plot(activity.psths(1).window(1) : activity.psths(1).window(2), mean(data.(st), 1), "LineWidth", 4, "Color", "k");
    title(st);
    axis square; axis tight;
    xlim([-5, 30]);
    ylim([-0.25, 0.5]); ylabel("\DeltaF/F");

end
%% plot single-cell summary for 2D stimuli

% [~, idx] = sort(vertcat(spikes.kernels.peak), "descend");
% id = idx(7);
[~, idx] = sort(vertcat(spikes.stats.SNR), "descend");
id = idx(7);

f = cellplot2D(dFF, spikes, id);

f.Position = [122.1429 535.8571 1.9531e+03 528.5714];
% f.Position = [-158, 1101, 1715, 503];

%% plot population summary for 2D stimuli

f = popplot2D(dFF, spikes);

f.Position = [1.1016e+03 466.7143 741.1429 647.4286];

%% plot single-cell summary for 3D stimuli

[~, idx] = sort(vertcat(spikes.kernels.peak), "descend");
id = 116; % idx(2);

f = figure; set(gcf, "Color", "w"); tiledlayout(2, 7, "TileSpacing", "tight");

%%%% dFF
% % trace
% nexttile([1, 2]); plot(dFF.cells(id).data, "LineWidth", 0.5, "Color", "k"); title("\DeltaF / F"); axis tight;
% PSTH
nexttile; hold on; plot_stim_resp(dFF.cells(id).data, frame_on, dFF.viz.window);
% contrast/size kernel
nexttile; imagesc(squeeze(dFF.kernels(id).kernsmooth(:, dFF.kernels(id).peakidx(2), :)));
xticks(dFF.viz.ticks.disc{3}); xticklabels(dFF.viz.ticklabels{3}); % xlabel(dFF.viz.labels{3});
yticks(dFF.viz.ticks.disc{1}); yticklabels(dFF.viz.ticklabels{1}); ylabel(dFF.viz.labels{1});
% contrast/ori kernel
nexttile; imagesc(squeeze(dFF.kernels(id).kernsmooth(:, :, dFF.kernels(id).peakidx(3))));
xticks(dFF.viz.ticks.disc{2}); xticklabels(dFF.viz.ticklabels{2}); % xlabel(dFF.viz.labels{2});
yticks(dFF.viz.ticks.disc{1}); yticklabels(dFF.viz.ticklabels{1}); ylabel(dFF.viz.labels{1});
% size/ori kernel
nexttile; imagesc(squeeze(dFF.kernels(id).kernsmooth(dFF.kernels(id).peakidx(1), :, :))'); colorbar;
xticks(dFF.viz.ticks.disc{2}); xticklabels(dFF.viz.ticklabels{2}); % xlabel(dFF.viz.labels{2});
yticks(dFF.viz.ticks.disc{3}); yticklabels(dFF.viz.ticklabels{3}); ylabel(dFF.viz.labels{3});
% contrast curve
nexttile; plot(dFF.viz.ticks.cont{1}, dFF.kernels(id).curve{1}, "-ok"); axis tight;
xticks(dFF.viz.ticks.cont{1}); xticklabels(dFF.viz.ticklabels{1}); xlabel(dFF.viz.labels{1});
% ori curve
nexttile; plot(dFF.viz.ticks.cont{2}, dFF.kernels(id).curve{2}, "-ok"); axis tight;
xticks(dFF.viz.ticks.cont{2}); xticklabels(dFF.viz.ticklabels{2}); xlabel(dFF.viz.labels{2});
% size curve
nexttile; plot(dFF.viz.ticks.cont{3}, dFF.kernels(id).curve{3}, "-ok"); axis tight;
xticks(dFF.viz.ticks.cont{3}); xticklabels(dFF.viz.ticklabels{3}); xlabel(dFF.viz.labels{3});

%%%% spikes
% % trace
% nexttile([1, 2]); plot(spikes.cells(id).data, "LineWidth", 0.5, "Color", "k"); title("\DeltaF / F"); axis tight;
% PSTH
nexttile; hold on; plot_stim_resp(spikes.cells(id).data, frame_on, spikes.viz.window);
% contrast/size kernel
nexttile; imagesc(squeeze(spikes.kernels(id).kernsmooth(:, spikes.kernels(id).peakidx(2), :))); pbaspect([flip(size(dFF.kernels(id).kern, [1, 3])), 1]);
xticks(spikes.viz.ticks.disc{3}); xticklabels(spikes.viz.ticklabels{3}); xlabel(spikes.viz.labels{3});
yticks(spikes.viz.ticks.disc{1}); yticklabels(spikes.viz.ticklabels{1}); ylabel(spikes.viz.labels{1});
% contrast/ori kernel
nexttile; imagesc(squeeze(spikes.kernels(id).kernsmooth(:, :, spikes.kernels(id).peakidx(3)))); pbaspect([flip(size(dFF.kernels(id).kern, [1, 2])), 1]);
xticks(spikes.viz.ticks.disc{2}); xticklabels(spikes.viz.ticklabels{2}); xlabel(spikes.viz.labels{2});
yticks(spikes.viz.ticks.disc{1}); yticklabels(spikes.viz.ticklabels{1}); ylabel(spikes.viz.labels{1});
% size/ori kernel
nexttile; imagesc(squeeze(spikes.kernels(id).kernsmooth(spikes.kernels(id).peakidx(1), :, :))'); pbaspect([flip(size(dFF.kernels(id).kern, [3, 2])), 1]); colorbar; 
xticks(spikes.viz.ticks.disc{2}); xticklabels(spikes.viz.ticklabels{2}); xlabel(spikes.viz.labels{2});
yticks(spikes.viz.ticks.disc{3}); yticklabels(spikes.viz.ticklabels{3}); ylabel(spikes.viz.labels{3});
% contrast curve
nexttile; plot(spikes.viz.ticks.cont{1}, spikes.kernels(id).curve{1}, "-ok"); axis tight;
xticks(spikes.viz.ticks.cont{1}); xticklabels(spikes.viz.ticklabels{1}); xlabel(spikes.viz.labels{1});
% ori curve
nexttile; plot(spikes.viz.ticks.cont{2}, spikes.kernels(id).curve{2}, "-ok"); axis tight;
xticks(spikes.viz.ticks.cont{2}); xticklabels(spikes.viz.ticklabels{2}); xlabel(spikes.viz.labels{2});
% size curve
nexttile; plot(spikes.viz.ticks.cont{3}, spikes.kernels(id).curve{3}, "-ok"); axis tight;
xticks(spikes.viz.ticks.cont{3}); xticklabels(spikes.viz.ticklabels{3}); xlabel(spikes.viz.labels{3});

f.Position = [1 589.5714 2.1914e+03 434.2857];
% f.Position = [-158, 1101, 1715, 503];

%% plot population summary for 3D stimuli

peakvals = vertcat(spikes.kernels.peakval);
peakvals = peakvals(peakvals(:, 1) > 0, :); % remove cells preferring zero contrast

figure; tiledlayout(1, 2, "TileSpacing", "tight");
nexttile; histogram(log10(peakvals(:, 1)), spikes.viz.binedges{1});
xticks(spikes.viz.ticks.cont{1}); xticklabels(spikes.viz.ticklabels{1}); xlabel(spikes.viz.labels{1});
xline(median(log10(peakvals(:, 1))), "--r", "LineWidth", 2);
setStyle();
ylabel("# of Cells");
% nexttile; histogram(peakvals(:, 2), spikes.viz.binedges{2});
% xticks(spikes.viz.xticks.cont); xticklabels(spikes.viz.xticklabels); xlabel(spikes.viz.xlabel);
nexttile; histogram(peakvals(:, 3), spikes.viz.binedges{3});
xticks(spikes.viz.ticks.cont{3}); xticklabels(spikes.viz.ticklabels{3}); xlabel(spikes.viz.labels{3});
xline(median(peakvals(:, 3)), "--r", "LineWidth", 2);
setStyle();

estvals = vertcat(spikes.kernels.estval);
estvals = estvals(estvals(:, 1) > 0, :); % remove cells preferring zero contrast

figure; tiledlayout(1, 2, "TileSpacing", "tight");
nexttile; histogram(log10(estvals(:, 1)), spikes.viz.binedges{1});
xticks(spikes.viz.ticks.cont{1}); xticklabels(spikes.viz.ticklabels{1}); xlabel(spikes.viz.labels{1});
xline(median(log10(estvals(:, 1))), "--r", "LineWidth", 2);
setStyle();
ylabel("# of Cells");
% nexttile; histogram(peakvals(:, 2), spikes.viz.binedges{2});
% xticks(spikes.viz.xticks.cont); xticklabels(spikes.viz.xticklabels); xlabel(spikes.viz.xlabel);
nexttile; histogram(estvals(:, 3), spikes.viz.binedges{3});
xticks(spikes.viz.ticks.cont{3}); xticklabels(spikes.viz.ticklabels{3}); xlabel(spikes.viz.labels{3});
xline(median(estvals(:, 3)), "--r", "LineWidth", 2);
setStyle();

%% plot for 4D stimuli

[~, idx] = sort(vertcat(spikes.kernels.peak), "descend");
id = idx(7);

f = figure; set(gcf, "Color", "w"); tiledlayout(2, 4);
% traces
nexttile([1, 2]); plot(dFF.cells(id).data, "LineWidth", 0.5, "Color", "k"); title("\DeltaF / F"); axis tight;
nexttile; hold on; plot_stim_resp(dFF.cells(id).data, frame_on, dFF.viz.window);
nexttile; imagesc(squeeze(mean(dFF.kernels(id).kern, [3, 4]))); colorbar;
xticks(dFF.viz.xticks.disc); xticklabels(dFF.viz.xticklabels); xlabel(dFF.viz.xlabel);
yticks(dFF.viz.yticks.disc); yticklabels(dFF.viz.yticklabels); ylabel(dFF.viz.ylabel);
% spikes
nexttile([1, 2]); plot(spikes.cells(id).data, "LineWidth", 0.5, "Color", "k"); title("Inferred Spikes"); axis tight;
nexttile; hold on; plot_stim_resp(spikes.cells(id).data, frame_on, spikes.viz.window);
nexttile; imagesc(squeeze(mean(spikes.kernels(id).kern, [3, 4]))); colorbar;
xticks(spikes.viz.xticks.disc); xticklabels(spikes.viz.xticklabels); xlabel(spikes.viz.xlabel);
yticks(spikes.viz.yticks.disc); yticklabels(spikes.viz.yticklabels); ylabel(spikes.viz.ylabel);

% f.Position = [769.5714 530.1429 1.3291e+03 576.5714];
