
%%%% add tools to path
addpath(genpath(fullfile("E:/Imaging/code/tools/")));

expsts = [];
path = "E:/Imaging/";
subjects = ["pv_vip_01", "pv_vip_04", "pv_vip_06", "pv_vip_07"]; % 
celltypes = ["Pyr", "Pyr", "Pyr", "Sst", "Sst", "Sst", "Sst", "Pyr"]; % 
recordings = {["6", "7"], ["1"], ["5", "6", "7", "8"], ["1"]}; % 
stims = ["trinoise", "randorisf", "battery1", "battery2", "battery3", "battery4"];
sid = 1:6;
stim_radius = 20; pix_per_deg = 14.3258;
eye_radius = 4;

for sj = 1:numel(subjects)

    for rc = 1:numel(recordings{sj})

        sessions = recordings{sj}(rc) + ["10", "20", "30", "40", "50", "60"];

        for st = sid
        
            tic;

            % load calcium data
            [rtdata, stdata, ops, stat] = load_calcium_data(path, subjects(sj), sessions(st), "000");
            
            % load stimulus data
            [frame_on, frame_off, params, stim_center] = load_stimulus_data(path, subjects(sj), sessions(st), "000");
            
            % load quadrature data
            quad = load_quad_data(path, subjects(sj), sessions(st), "000");
            
            % load eye data
            [eye, imgs] = load_eye_data(path, subjects(sj), sessions(st), "000");
            if stims(st) == "trinoise"
                eye_center = NaN; 
            else
                eye_center = expsts.(subjects(sj) + "_" + recordings{sj}(rc)).trinoise.spikes.eye.center;
            end
            
            % create experiment structure
            [expsts.(subjects(sj) + "_" + recordings{sj}(rc)).(stims(st)).spikes, ...
             expsts.(subjects(sj) + "_" + recordings{sj}(rc)).(stims(st)).dFF] = expstruct(stims(st), ...
                                                                                           rtdata, stdata, ops, stat, ...
                                                                                           frame_on, frame_off, params, ...
                                                                                           stim_center, pix_per_deg, ...
                                                                                           eye, imgs, eye_center, eye_radius, ...
                                                                                           quad);
        
            fprintf("Loaded %s for recording %s from %s...\n", stims(st), recordings{sj}(rc), subjects(sj));
            toc;
        
        end

    end

end

%% make population summary

cells = structfun(@(x) filter_cells(x), expsts, "UniformOutput", false);

f = figure; tiledlayout(4, 20);

colors = {[0.7, 0.7, 0.7], [0.15, 0.5, 0.8]};

% randorisf histograms incl. ori, sf, SNR (and thr), F1F0
nexttile(1, [2, 2]); hold on; feature_hist_cats(expsts, "randorisf", "spikes", celltypes, colors, "SNR_ROSF_spikes", cells, "SNR", true);
nexttile(41, [2, 2]); hold on; feature_hist_cats(expsts, "randorisf", "spikes", celltypes, colors, "F1F0_ROSF_spikes", cells, "F1F0", false);
nexttile(3, [2, 2]); hold on; feature_hist_cats(expsts, "randorisf", "spikes", celltypes, colors, "Ori_ROSF_spikes", cells, "Orientation", false);
nexttile(43, [2, 2]); hold on; feature_hist_cats(expsts, "randorisf", "spikes", celltypes, colors, "SF_ROSF_spikes", cells, "Spatial Frequency", false);

% battery1 histograms, incl. dir, sf, DSI
nexttile(5, [2, 2]); hold on; feature_hist_cats(expsts, "battery1", "spikes", celltypes, colors, "Dir_B1_spikes", cells, "Direction", false);
nexttile(45, [2, 2]); hold on; feature_hist_cats(expsts, "battery1", "spikes", celltypes, colors, "SF_B1_spikes", cells, "Spatial Frequency", false);
nexttile(7, [2, 2]); hold on; feature_hist_cats(expsts, "battery1", "spikes", celltypes, colors, "DSI_B1_spikes", cells, "DSI", false);

% battery2 histograms, incl. tf
nexttile(47, [2, 2]); hold on; feature_hist_cats(expsts, "battery2", "spikes", celltypes, colors, "TF_B2_spikes", cells, "Temporal Frequency", false);

% battery3 histograms, incl. contrast, size
nexttile(9, [2, 2]); hold on; feature_hist_cats(expsts, "battery3", "spikes", celltypes, colors, "Ctrst_B3_spikes", cells, "Contrast", false);
nexttile(49, [2, 2]); hold on; feature_hist_cats(expsts, "battery3", "spikes", celltypes, colors, "Size_B3_spikes", cells, "Size", false);
% running/stationary size tuning curve
nexttile(11, [2, 2]); axis square; axis tight; hold on; size_running_cats(expsts, "battery3", "dFF", celltypes, colors, cells);

% cross vs. iso
nexttile(51, [2, 2]); hold on; ctxmod_scatter_cats(expsts, "battery4", "dFF", ["Cross", "Iso_All"], celltypes, colors, cells, 3); % use 20 degrees

% cross/iso tuning curves at different sizes
nexttile(13, [2, 2]); axis square; axis tight; hold on; ctxmod_curve_cats(expsts, "battery4", "dFF", ["Iso", "Cross"], celltypes, colors, cells);

% center/surround tuning curves at different sizes
nexttile(53, [2, 2]); axis square; axis tight; hold on; ctxmod_curve_cats(expsts, "battery4", "dFF", ["Surround", "Center"], celltypes, colors, cells);

% battery4 histograms, incl. modulation index, ITI
nexttile(15, [2, 2]); hold on; feature_hist_cats(expsts, "battery4", "dFF", celltypes, colors, "cri_mod_mean_B4_dFF", cells, "Cross/Iso Modulation Index", false);
nexttile(55, [2, 2]); hold on; feature_hist_cats(expsts, "battery4", "dFF", celltypes, colors, "ci_mod_mean_B4_dFF", cells, "Center/Iso Modulation Index", false);
nexttile(17, [2, 2]); hold on; feature_hist_cats(expsts, "battery4", "dFF", celltypes, colors, "cs_size_corr_B4_dFF", cells, "C/S Size Tuning Correlation", false);
% contrast-response function
nexttile(57, [2, 2]); axis square; axis tight; hold on;
contrast_response(expsts, "battery3", "dFF", celltypes, colors, cells);

set(f, "Color", "w");
f.Position = [41.5714 434.1429 2.1103e+03 639.4285];

%% cluster and classify cells

[tbl, ctypes, ids] = build_table(expsts, celltypes, cells);
clrs = vertcat(colors{:});
clrs = clrs((ctypes == "Sst") + 1, :);

vars = ["cri_mod_mean_B4_dFF", "ci_mod_mean_B4_dFF", "cs_size_corr_B4_dFF", ...
        "Ori_ROSF_spikes", "SF_ROSF_spikes", "F1F0_ROSF_spikes", "CV_ROSF_spikes", ...
        "SF_B1_spikes", "Dir_B1_spikes", "DSI_B1_spikes", ...
        "TF_B2_spikes", ...
        "Ctrst_B3_spikes", "Size_B3_spikes"];

% vars = ["cri_mod_mean_B4_dFF", "ci_mod_mean_B4_dFF", "cs_size_corr_B4_dFF"];

tblsub = tbl(:, vars);
arr = table2array(tblsub);
% remove NaN
nanidx = any(isnan(arr), 2);
arr(nanidx, :) = [];
clrs(nanidx, :) = [];
tblz = zscore(arr, 0, 1);
X = arr;
[~, S] = pca(tblz);

f = figure; tiledlayout(1, 3);
nexttile; scatter(S(:, 1), S(:, 2), 20, clrs, "filled"); axis square; xlabel("PC 1"); ylabel("PC 2");
Y = tsne(S);
nexttile; hold on;
pyr_idx = ctypes == "Pyr"; sst_idx = ctypes == "Sst";
scatter(Y(pyr_idx(~nanidx), 1), Y(pyr_idx(~nanidx), 2), 20, "filled", "MarkerFaceColor", colors{1}, "MarkerEdgeColor", colors{1});
scatter(Y(sst_idx(~nanidx), 1), Y(sst_idx(~nanidx), 2), 20, "filled", "MarkerFaceColor", colors{2}, "MarkerEdgeColor", colors{2});
axis square; xlabel("t-SNE 1"); ylabel("t-SNE 2");
legend("Pyr", "SST", "Box", "off");

% classifier (balanced classes)

pyr_recs = ["pv_vip_01_6", "pv_vip_01_7", "pv_vip_04_1", "pv_vip_07_1"];
sst_recs = ["pv_vip_06_5", "pv_vip_06_6", "pv_vip_06_7", "pv_vip_06_8"];

perf = cell(numel(pyr_recs) * numel(sst_recs), 2);
err = nan(numel(pyr_recs) * numel(sst_recs), 2);

cv_idx = 1;
for sim = 1:10
    for pyr = pyr_recs
        for sst = sst_recs
            X = tblz;
            y = ctypes(~nanidx) == "Sst";
            dset = ids(~nanidx, 1);
            all_idx = 1:size(X, 1);
            train_idx_pyr = ismember(dset, setdiff(pyr_recs, pyr));
            train_idx_sst = ismember(dset, setdiff(sst_recs, sst));
            test_idx_pyr = ismember(dset, pyr); if nnz(test_idx_pyr) == 0, break; end
            test_idx_sst = ismember(dset, sst); if nnz(test_idx_sst) == 0, continue; end
            [train_idx_pyr, train_idx_sst] = balance_classes(train_idx_pyr, train_idx_sst);
            [test_idx_pyr, test_idx_sst] = balance_classes(test_idx_pyr, test_idx_sst);
            train_idx = train_idx_pyr | train_idx_sst;
            test_idx = test_idx_pyr | test_idx_sst;
            mdl = fitcdiscr(X(train_idx, :), y(train_idx));
            y_pred = predict(mdl, X(test_idx, :));
            X_perm = X(test_idx, :);
            y_perm = predict(mdl, X_perm(randperm(size(X_perm, 1)), :));
            perf{cv_idx, 1} = classperf(y(test_idx), y_perm);
            perf{cv_idx, 2} = classperf(y(test_idx), y_pred);
            err(cv_idx, :) = [perf{cv_idx, 1}.CorrectRate, perf{cv_idx, 2}.CorrectRate];
            cv_idx = cv_idx + 1;
        end
    end
end

nexttile; hold on; axis square;
histogram(err(:, 1), 0:0.1:1, "Normalization", "probability", "FaceColor", "white", "EdgeColor", "black");
histogram(err(:, 2), 0:0.1:1, "Normalization", "probability", "FaceColor", [0.7, 0.7, 0.7], "EdgeColor", [0.7, 0.7, 0.7]);
xlabel("Cross-Validated Accuracy"); ylabel("Fraction of Folds");
legend("Shuffled", "Observed", "box", "off", "Location", "northwest");

set(f, "Color", "w");
f.Position = [645 507.8571 718.2857 306.2857];

%% check on individual cells

expst = expsts.pv_vip_06_5;
id = 15;

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

%% functions

function cells = filter_cells(expsts)
    % ON domain criteria
    on_dist = vertcat(expsts.trinoise.spikes.stats.on_dist_fit) < 10;
    on_snr = vertcat(expsts.trinoise.spikes.stats.snr_on) > 5;
    on_ve = vertcat(expsts.trinoise.spikes.stats.cc_on) > 0.5;
    on_centered = on_dist & on_snr & on_ve;
    on_noise = ~(on_snr & on_ve);
    % OFF domain criteria
    off_dist = vertcat(expsts.trinoise.spikes.stats.off_dist_fit) < 10;
    off_snr = vertcat(expsts.trinoise.spikes.stats.snr_off) > 5;
    off_ve = vertcat(expsts.trinoise.spikes.stats.cc_off) > 0.5;
    off_centered = off_dist & off_snr & off_ve;
    off_noise = ~(off_snr & off_ve);
    % visually responsive
    SNR = vertcat(expsts.randorisf.spikes.stats.SNR) > expsts.randorisf.spikes.calc.SNR_thr;
    % intersect
    cells = SNR & ((on_centered & off_centered) | (on_centered & off_noise) | (off_centered & on_noise));
end

function feature_hist_cats(expsts, session, dtype, celltypes, colors, featname, idx, t, lgd)
    recs = fieldnames(expsts);
    ctypes = unique(celltypes);
    data = cell2struct(cell(numel(ctypes),1), ctypes);
    for r = 1:numel(recs)
        tbl = expsts.(recs{r}).(session).(dtype).table;
        data.(celltypes(r)) = vertcat(data.(celltypes(r)), tbl.(featname)(idx.(recs{r})));
    end
    feat = ismember(tbl.Properties.VariableNames, featname);
    scale = string(tbl.Properties.CustomProperties.FeatureScale{feat});
    binedges = tbl.Properties.CustomProperties.BinEdges{feat};
    ticks = tbl.Properties.CustomProperties.Ticks{feat};
    ticklabels = tbl.Properties.CustomProperties.TickLabels{feat};
    for c = 1:numel(ctypes)
        if scale == "linear" || scale == "polar"
            histogram(data.(ctypes(c)), binedges, "Normalization", "probability", "FaceColor", colors{c}, "EdgeColor", colors{c});
        elseif scale == "log"
            histogram(log10(data.(ctypes(c))), binedges, "Normalization", "probability", "FaceColor", colors{c}, "EdgeColor", colors{c});
        end
    end
    xticks(ticks); xticklabels(ticklabels);
    title(t, "Interpreter", "none");
    axis square; axis tight;
    if lgd == true
        legend(ctypes, "Box", "off", "Location", "northeast");
    end
end

function size_running_cats(expsts, session, dtype, celltypes, colors, idx)
    recs = fieldnames(expsts);
    ctypes = unique(celltypes);
    data = cell2struct(repmat({struct("stationary", [], "running", [])}, [numel(ctypes), 1]), ctypes);
    for r = 1:numel(recs)
        stationary = vertcat(expsts.(recs{r}).(session).(dtype).stats.size_stationary);
        running = vertcat(expsts.(recs{r}).(session).(dtype).stats.size_running);
        data.(celltypes(r)).stationary = vertcat(data.(celltypes(r)).stationary, stationary(idx.(recs{r}), :));
        data.(celltypes(r)).running = vertcat(data.(celltypes(r)).running, running(idx.(recs{r}), :));
    end
    % plot
    for c = 1:numel(ctypes)
        maxval = max(horzcat(mean(data.(ctypes(c)).stationary, 1), mean(data.(ctypes(c)).running, 1)), [], "all");
        snorm = mean(data.(ctypes(c)).stationary, 1, "omitnan") ./ maxval;
        rnorm = mean(data.(ctypes(c)).running, 1, "omitnan") ./ maxval;
        plot(expsts.(recs{1}).(session).(dtype).viz.ticks.cont{3}, snorm, "Color", colors{c});
        scatter(expsts.(recs{1}).(session).(dtype).viz.ticks.cont{3}, snorm, "filled", "MarkerFaceColor", [1, 1, 1], "MarkerEdgeColor", colors{c});
        plot(expsts.(recs{1}).(session).(dtype).viz.ticks.cont{3}, rnorm, "Color", colors{c});
        scatter(expsts.(recs{1}).(session).(dtype).viz.ticks.cont{3}, rnorm, "filled", "MarkerFaceColor", colors{c}, "MarkerEdgeColor", colors{c});
    end
    xticks(expsts.(recs{1}).(session).(dtype).viz.ticks.cont{3}); xticklabels(expsts.(recs{1}).(session).(dtype).viz.ticklabels{3}); xlabel(expsts.(recs{1}).(session).(dtype).viz.labels{3});
    ylabel("\DeltaF/F (Normalized)");
    legend("", "Stationary", "", "Running", "", "", "", "", "Box", "off", "Location", "northoutside", "Orientation", "horizontal");
end

function ctxmod_curve_cats(expsts, session, dtype, cats, celltypes, colors, idx)
    recs = fieldnames(expsts);
    ctypes = unique(celltypes);
    data = cell2struct(repmat({struct(cats{1}, [], cats{2}, [])}, [numel(ctypes), 1]), ctypes);
    for r = 1:numel(recs)
        cat1 = vertcat(expsts.(recs{r}).(session).(dtype).stats.(lower(cats{1}) + "_mean"));
        cat2 = vertcat(expsts.(recs{r}).(session).(dtype).stats.(lower(cats{2}) + "_mean"));
        data.(celltypes(r)).(cats{1}) = vertcat(data.(celltypes(r)).(cats{1}), cat1(idx.(recs{r}), :));
        data.(celltypes(r)).(cats{2}) = vertcat(data.(celltypes(r)).(cats{2}), cat2(idx.(recs{r}), :));
    end
    % plot
    for c = 1:numel(ctypes)
        maxval = max(horzcat(mean(data.(ctypes(c)).(cats{1}), 1), mean(data.(ctypes(c)).(cats{2}), 1)), [], "all");
        norm1 = mean(data.(ctypes(c)).(cats{1}), 1, "omitnan") ./ maxval;
        norm2 = mean(data.(ctypes(c)).(cats{2}), 1, "omitnan") ./ maxval;
        plot(expsts.(recs{1}).(session).(dtype).viz.ticks.cont{5}, norm1, "Color", colors{c});
        scatter(expsts.(recs{1}).(session).(dtype).viz.ticks.cont{5}, norm1, "filled", "MarkerFaceColor", [1, 1, 1], "MarkerEdgeColor", colors{c});
        plot(expsts.(recs{1}).(session).(dtype).viz.ticks.cont{5}, norm2, "Color", colors{c});
        scatter(expsts.(recs{1}).(session).(dtype).viz.ticks.cont{5}, norm2, "filled", "MarkerFaceColor", colors{c}, "MarkerEdgeColor", colors{c});
    end
    xticks(expsts.(recs{1}).(session).(dtype).viz.ticks.cont{5}); xticklabels(expsts.(recs{1}).(session).(dtype).viz.ticklabels{5}); xlabel(expsts.(recs{1}).(session).(dtype).viz.labels{5});
    ylabel("\DeltaF/F (Normalized)");
    legend("", cats{1}, "", cats{2}, "", "", "", "", "Box", "off", "Location", "northoutside", "Orientation", "horizontal");
end

function ctxmod_scatter_cats(expsts, session, dtype, cats, celltypes, colors, idx, sz_idx)
    recs = fieldnames(expsts);
    ctypes = unique(celltypes);
    data = cell2struct(cell(numel(ctypes),1), ctypes);
    for r = 1:numel(recs)
        cat1 = vertcat(expsts.(recs{r}).(session).(dtype).stats.(lower(cats{1}) + "_mean"));
        cat2 = vertcat(expsts.(recs{r}).(session).(dtype).stats.(lower(cats{2}) + "_mean"));
        catcat = horzcat(cat1(:, sz_idx), cat2(:, sz_idx));
        data.(celltypes(r)) = vertcat(data.(celltypes(r)), catcat(idx.(recs{r}), :));
    end
    % scatter
    for c = 1:numel(ctypes)
        X = data.(ctypes{c})(:, 1);
        Y = data.(ctypes{c})(:, 2);
        xline(0, "--k", "LineWidth", 1);
        yline(0, "--k", "LineWidth", 1);
        plot(-1:1, -1:1, "--k", "LineWidth", 1);
        X(X > 1) = 1;
        Y(Y > 1) = 1;
        scatter(X, Y, 10, colors{c}, "filled", "MarkerFaceAlpha", 0.75, "MarkerEdgeAlpha", 0.75);
    end
    xlim([-0.1, 0.6]); xticks(-0.1:0.1:0.6);
    ylim([-0.1, 0.6]); yticks(-0.1:0.1:0.6);
    xlabel([erase(cats{1}, "_All"), ' [\DeltaF/F]']); ylabel([erase(cats{2}, "_All"), ' [\DeltaF/F]']);
    axis square;
end

function [tbl, ctypes, ids] = build_table(expsts, celltypes, idx)
    recs = fieldnames(expsts);
    tbl = table();
    ctypes = [];
    ids = [];
    for r = 1:numel(recs)
        sessions = setdiff(fieldnames(expsts.(recs{r})), "trinoise");
        tblrec = table();
        for s = 1:numel(sessions)
            tblrec = horzcat(tblrec, horzcat(expsts.(recs{r}).(sessions{s}).spikes.table, ...
                                             expsts.(recs{r}).(sessions{s}).dFF.table));
        end
        tbl = vertcat(tbl, tblrec(idx.(recs{r}), :));
        ctypes = vertcat(ctypes, repmat(celltypes(r), [nnz(idx.(recs{r})), 1]));
        ids = vertcat(ids, horzcat(repmat(recs(r), [nnz(idx.(recs{r})), 1]), num2cell(find(idx.(recs{r})))));
    end
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

function contrast_response(expsts, session, dtype, celltypes, colors, idx)
    recs = fieldnames(expsts);
    ctypes = unique(celltypes);
    data = cell2struct(cell(numel(ctypes),1), ctypes);
    for r = 1:numel(recs)
        cat = vertcat(expsts.(recs{r}).(session).(dtype).kernels.curve);
        cat = vertcat(cat{:, 1});
        data.(celltypes(r)) = vertcat(data.(celltypes(r)), cat(idx.(recs{r}), :));
    end
    % plot lines
    for c = 1:numel(ctypes)
        % maxval = max(horzcat(mean(data.(ctypes(c)).(cats{1}), 1), mean(data.(ctypes(c)).(cats{2}), 1)), [], "all");
        X = log10(expsts.(recs{1}).(session).(dtype).dimvals{1});
        % X(1) = X(2) - mean(diff(X(2:end)));
        M = mean(data.(ctypes(c)), 1, "omitnan");
        Y = log10(M - min(M) + 0.01);
        plot(X(2:end), Y(2:end), "Color", colors{c});
        scatter(X(2:end), Y(2:end), "filled", "MarkerFaceColor", colors{c}, "MarkerEdgeColor", colors{c});
    end
    % xticks(X); % xticklabels(expsts.(recs{1}).(session).(dtype).dimvals{1});
    % ylim([-0.1, 0.6]); yticks(-0.1:0.1:0.6);
    xlabel('log(Contrast)'); ylabel('log(\DeltaF/F)');
    % axis square;
end

function [idx_class_1, idx_class_2] = balance_classes(idx_class_1, idx_class_2)
    if nnz(idx_class_1) > nnz(idx_class_2)
        idx_1 = find(idx_class_1);
        rand_idx_1 = idx_1(randperm(numel(idx_1)));
        idx_class_1(rand_idx_1(1:nnz(idx_class_1) - nnz(idx_class_2))) = 0;
    elseif nnz(idx_class_1) < nnz(idx_class_2)
        idx_2 = find(idx_class_2);
        rand_idx_2 = idx_2(randperm(numel(idx_2)));
        idx_class_2(rand_idx_2(1:nnz(idx_class_2) - nnz(idx_class_1))) = 0;
    end
end
