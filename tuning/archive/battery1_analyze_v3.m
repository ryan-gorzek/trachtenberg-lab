
% add tools to path
addpath(genpath(fullfile("E:/Imaging/code/tools/")));

% load calcium data (from realtime for now)
[rtdata, stdata] = load_calcium_data("pv_vip_00", "340", "001");
% remove bad cells
rtdata(:,isnan(mean(rtdata))) = [];
stdata(:,isnan(mean(stdata))) = [];

% load stimulus data (2D for now)
[frame_on, frame_off, params] = load_stimulus_data("pv_vip_00", "340", "001");

% % load quadrature data
% quad = load_quadrature_data("pv_vip_00", "310", "000");
% 
% % load eye data
% eye = load_eye_data("pv_vip_00", "310", "000");

%%%% Cells Structure
%%%% allowing for streamlined access, sorting, plotting
% traces (various forms), spikes
% kernels
% max/min kernels, optimal values
% correlation with running/eye
% signal, noise, variance, etc.

% create experiment structure
[spikes, dFF, ~] = expstruct("battery4", rtdata, stdata, frame_on, frame_off, params);

%% plot for 2D stimuli

[~, idx] = sort(vertcat(spikes.kernels.peak), "descend");
id = idx(10);

f = figure; set(gcf, "Color", "w"); tiledlayout(2, 4);
% traces
nexttile([1, 2]); plot(dFF.cells(id).data, "LineWidth", 0.5, "Color", "k"); title("\DeltaF / F"); axis tight;
nexttile; hold on; plot_stim_resp(dFF.cells(id).data, frame_on, dFF.viz.window);
nexttile; imagesc(dFF.kernels(id).kern); colorbar;
xticks(dFF.viz.xticks.disc); xticklabels(dFF.viz.xticklabels); xlabel(dFF.viz.xlabel);
yticks(dFF.viz.yticks.disc); yticklabels(dFF.viz.yticklabels); ylabel(dFF.viz.ylabel);
% spikes
nexttile([1, 2]); plot(spikes.cells(id).data, "LineWidth", 0.5, "Color", "k"); title("Inferred Spikes"); axis tight;
nexttile; hold on; plot_stim_resp(spikes.cells(id).data, frame_on, spikes.viz.window);
nexttile; imagesc(spikes.kernels(id).kern); colorbar;
xticks(spikes.viz.xticks.disc); xticklabels(spikes.viz.xticklabels); xlabel(spikes.viz.xlabel);
yticks(spikes.viz.yticks.disc); yticklabels(spikes.viz.yticklabels); ylabel(spikes.viz.ylabel);

f.Position = [769.5714 530.1429 1.3291e+03 576.5714];

%% plot for 3D stimuli

[~, idx] = sort(vertcat(spikes.kernels.peak), "descend");
id = idx(1);

f = figure; set(gcf, "Color", "w"); tiledlayout(2, 4);
% traces
nexttile([1, 2]); plot(dFF.cells(id).data, "LineWidth", 0.5, "Color", "k"); title("\DeltaF / F"); axis tight;
nexttile; hold on; plot_stim_resp(dFF.cells(id).data, frame_on, dFF.viz.window);
nexttile; imagesc(squeeze(mean(dFF.kernels(id).kern, 2))); colorbar;
xticks(dFF.viz.zticks.disc); xticklabels(dFF.viz.zticklabels); xlabel(dFF.viz.zlabel);
yticks(dFF.viz.yticks.disc); yticklabels(dFF.viz.yticklabels); ylabel(dFF.viz.ylabel);
% spikes
nexttile([1, 2]); plot(spikes.cells(id).data, "LineWidth", 0.5, "Color", "k"); title("Inferred Spikes"); axis tight;
nexttile; hold on; plot_stim_resp(spikes.cells(id).data, frame_on, spikes.viz.window);
nexttile; imagesc(squeeze(mean(spikes.kernels(id).kern, 2))); colorbar;
xticks(spikes.viz.zticks.disc); xticklabels(spikes.viz.zticklabels); xlabel(spikes.viz.zlabel);
yticks(spikes.viz.yticks.disc); yticklabels(spikes.viz.yticklabels); ylabel(spikes.viz.ylabel);

f.Position = [769.5714 530.1429 1.3291e+03 576.5714];

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

f.Position = [769.5714 530.1429 1.3291e+03 576.5714];
