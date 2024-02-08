
%%%% add tools to path
addpath(genpath(fullfile("E:/Imaging/code/tools/")));
% addpath(genpath(fullfile("/Users/ryan.gorzek/Imaging/code/tools/")));

path = "E:/Imaging";
stimdict = dictionary(["01", "40", "50"], ["randorisf", "battery4", "battery5"]);
subjects = ["pv_vip_01"];

tbl = table(); SNR_thr = [];
for subj = subjects

    tbl_subj = table();
    f = string({dir(sprintf("%s/%s/", path, subj)).name});
    f_idx = strncmp(subj, f, numel(subj));
    n_idx = endsWith(f, "111"); % ignore the concatenated session
    stimf = f(f_idx & ~n_idx);

    for st = stimf

        stparts = strsplit(st, "_");
        sess = stparts(end - 1); run = stparts(end);
        stim = char(sess); stim = stimdict(stim(2:3));
        
        % load calcium data
        [rtdata, stdata] = load_calcium_data(path, subj, sess, run);
        
        % load stimulus data
        [frame_on, frame_off, params] = load_stimulus_data(path, subj, sess, run);
        
        % load quadrature data
        quad = load_quad_data(path, subj, sess, run);
        
        % load eye data
        eye = load_eye_data(path, subj, sess, run);
        
        % create experiment structure
        [spikes, dFF] = expstruct(stim, rtdata, stdata, frame_on, frame_off, params, eye, quad);
    
        if stim == "randorisf", SNR_thr = [SNR_thr, spikes.calc.SNR_thr]; end
    
        tbl_subj = horzcat(tbl_subj, dFF.table, spikes.table);
    
    end

    tbl_subj.subject = repmat(find(subjects == subj), [size(tbl_subj, 1), 1]);
    tbl_subj = tbl_subj(tbl_subj.SNR_ROSF_spikes > SNR_thr(subjects == subj), :); %%%% subset to remove non-responsive cells
    tbl = vertcat(tbl, tbl_subj);

end

%% get the table into its final form

%%%% add direction selectivity index

keepvars = ["subject", ...
            "cri_mod_max_B4_dFF", "cri_mod_mean_B4_dFF", "cs_mod_max_B4_dFF", "cs_mod_mean_B4_dFF", "ci_mod_max_B4_dFF", "ci_mod_mean_B4_dFF", ...
            "cri_mod_max_B5_dFF", "cri_mod_mean_B5_dFF", "cs_mod_max_B5_dFF", "cs_mod_mean_B5_dFF", "ci_mod_max_B5_dFF", "ci_mod_mean_B5_dFF", "ITI_B5_dFF"];
%             "SF_B1_spikes", "Dir_B1_spikes", "DSI_B1_spikes", ... % "eye_corr_B1_dFF", "quad_corr_B1_dFF", ...
%             "TF_B2_spikes", "Dir_B2_spikes", "DSI_B2_spikes", ... % "eye_corr_B2_dFF", "quad_corr_B2_dFF", ...
%             "Ctrst_B3_spikes", "Dir_B3_spikes", "Size_B3_spikes", "DSI_B3_spikes"]; % , ... % "eye_corr_B3_dFF", "quad_corr_B3_dFF", ...        

tbl_sub = tbl(:, keepvars);
% tbl_sub = tbl_sub_var(tbl_sub_var.SNR_ROSF_spikes > mean(SNR_thr), :);
% tbl_sub = tbl_sub_var;
tbl_sub(any(isnan(table2array(tbl_sub)), 2), :) = [];
colors = tbl_sub.subject;
tbl_sub_final = tbl_sub(:, keepvars(2:end));

%% examine all feature distributions

featnames = tbl_sub_final.Properties.VariableNames;
nfeats = numel(featnames);

% pooled
figure; tiledlayout(5, 3, "TileSpacing", "compact");
for feat = 1:nfeats
    scale = string(tbl_sub_final.Properties.CustomProperties.FeatureScale{feat});
    binedges = tbl_sub_final.Properties.CustomProperties.BinEdges{feat};
    ticks = tbl_sub_final.Properties.CustomProperties.Ticks{feat};
    ticklabels = tbl_sub_final.Properties.CustomProperties.TickLabels{feat};
    nexttile; hold on; setStyle(fontSize=10);
    if scale == "linear" || scale == "polar"
        histogram(tbl_sub_final.(featnames{feat}), binedges, "Normalization", "probability");
    elseif scale == "log"
        histogram(log10(tbl_sub_final.(featnames{feat})), binedges, "Normalization", "probability");
    end
    xticks(ticks); xticklabels(ticklabels);
    title(featnames{feat}, "Interpreter", "none");
    axis square; axis tight;
end

% separate
featnames = tbl_sub_final.Properties.VariableNames;
nfeats = numel(featnames);

% pooled
figure; tiledlayout(5, 3, "TileSpacing", "compact");
for feat = 1:nfeats
    scale = string(tbl_sub_final.Properties.CustomProperties.FeatureScale{feat});
    binedges = tbl_sub_final.Properties.CustomProperties.BinEdges{feat};
    ticks = tbl_sub_final.Properties.CustomProperties.Ticks{feat};
    ticklabels = tbl_sub_final.Properties.CustomProperties.TickLabels{feat};
    nexttile; hold on; setStyle(fontSize=10);
    for subj = 1:numel(unique(colors))
        idx = colors == subj;
        if scale == "linear" || scale == "polar"
            histogram(tbl_sub_final.(featnames{feat})(idx), binedges, "Normalization", "probability");
        elseif scale == "log"
            histogram(log10(tbl_sub_final.(featnames{feat})(idx)), binedges, "Normalization", "probability");
        end
    end
    xticks(ticks); xticklabels(ticklabels);
    title(featnames{feat}, "Interpreter", "none");
    axis square; axis tight;
end

%% cluster and plot cells in tsne space

tblz = zscore(table2array(tbl_sub_final), 0, 1);

[idx, C, sumdist] = kmeans(tblz, 2, 'Distance', 'cosine', ...
                                    'Display', 'final');

figure;
[silh, h] = silhouette(tblz, idx, 'cosine');
xlabel('Silhouette Value');
ylabel('Cluster');

Y = tsne(tblz);
figure; scatter(Y(:, 1), Y(:, 2), 30, colors);

%% examine correlations between variables

cmap = colorMap({[0, 0, 1], [1, 1, 1], [1, 0, 0]});
figure;
imagesc(corr(tblz)); colorbar; colormap(cmap); clim([-1, 1]);
set(gca, "TickLabelInterpreter", "None");
xticks(1:size(tblz, 2)); xticklabels(tbl_sub_final.Properties.VariableNames);
yticks(1:size(tblz, 2)); yticklabels(tbl_sub_final.Properties.VariableNames);

%% run pca

[coeff,score,latent,tsquared,explained,mu] = pca(tblz);

figure; scatter(score(:, 1), score(:, 2), 30, colors);

%% examine pupil and quadrature correlation distributions

idx_subj1 = tbl.subject == 1;
idx_subj2 = tbl.subject == 2;
figure; tiledlayout(1, 5);
binedges = -0.4:0.05:0.6;
nexttile; hold on; histogram(tbl.eye_corr_ROSF_dFF(idx_subj1), binedges); histogram(tbl.eye_corr_ROSF_dFF(idx_subj2), binedges);
nexttile; hold on; histogram(tbl.eye_corr_B1_dFF(idx_subj1), binedges); histogram(tbl.eye_corr_B1_dFF(idx_subj2), binedges);
nexttile; hold on; histogram(tbl.eye_corr_B2_dFF(idx_subj1), binedges); histogram(tbl.eye_corr_B2_dFF(idx_subj2), binedges);
nexttile; hold on; histogram(tbl.eye_corr_B3_dFF(idx_subj1), binedges); histogram(tbl.eye_corr_B3_dFF(idx_subj2), binedges);
nexttile; hold on; histogram(tbl.eye_corr_B4_dFF(idx_subj1), binedges); histogram(tbl.eye_corr_B4_dFF(idx_subj2), binedges);

figure; tiledlayout(1, 5);
binedges = -0.4:0.05:0.6;
nexttile; hold on; histogram(tbl.quad_corr_ROSF_dFF(idx_subj1), binedges); histogram(tbl.quad_corr_ROSF_dFF(idx_subj2), binedges);
nexttile; hold on; histogram(tbl.quad_corr_B1_dFF(idx_subj1), binedges); histogram(tbl.quad_corr_B1_dFF(idx_subj2), binedges);
nexttile; hold on; histogram(tbl.quad_corr_B2_dFF(idx_subj1), binedges); histogram(tbl.quad_corr_B2_dFF(idx_subj2), binedges);
nexttile; hold on; histogram(tbl.quad_corr_B3_dFF(idx_subj1), binedges); histogram(tbl.quad_corr_B3_dFF(idx_subj2), binedges);
nexttile; hold on; histogram(tbl.quad_corr_B4_dFF(idx_subj1), binedges); histogram(tbl.quad_corr_B4_dFF(idx_subj2), binedges);
