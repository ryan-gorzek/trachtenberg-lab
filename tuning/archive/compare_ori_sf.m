
%%%% Analyze orisf

subj_id = "pv_vip_00";
recordings = ["000_001"]; % , "210_000", "210_001"];

for r = 1:numel(recordings)
    exp_id = char(strcat(subj_id, "_", recordings(r)));
    fprintf("Processing %s...\n", exp_id);
    exp_path = char(strcat(exp_id, filesep));
    plane_path = strcat(exp_path, "suite2p", filesep, "plane0", filesep);
    % Generate .align, .segment, and .signals1 files. 
    sbxsuite2sbx(strcat(plane_path, "Fall"), strcat(exp_path, exp_id));
    % Generate .signals file.
    sbxf2spks(strcat(exp_path, exp_id));
    % Run orisf analysis.
    sbxorisf_new(strcat(exp_path, exp_id));
    % Load .orisf file and create table from structure array.
    stat_struct = load(strcat(exp_path, exp_id, ".orisf"), "-mat", "stat").stat;
end

%%

sper = 1920./(1.31.^(0:11));   % spatial period in pixels %1.395 gets to .32 cyc/deg
sper_pix = sper;
sper = sper / 15.25;           % 15.25 pixels per degree 
sf = 1./sper;                  % cycles/deg

ori = 0:20:170;

%%

response = vertcat(stat_struct.response);
noise = vertcat(stat_struct.noise);
noise_std = vertcat(stat_struct.noise_std);

figure; hold on;
plot(response');
plot(noise');
plot(noise' + 1 .* noise_std');

figure; hold on;
plot(response' - (noise' + 1 .* noise_std'), "-o", "LineWidth", 2);
yline(0, "--r", "LineWidth", 2);

idx_noise = response - (noise + 1 .* noise_std) <= 0;

%% plot the kernels

kernels = cat(3, stat_struct.kern);
kernels = kernels(:, :, idx_noise);
peak_ori = vertcat(stat_struct.ori_est);
peak_sf = vertcat(stat_struct.sf_est);

for kern = 1:size(kernels, 3)

    figure;
    imagesc(log10(sf), ori, kernels(:, :, kern));
    xlabel('Spatial Frequency (cycles/deg)');
    ylabel('Orientation (deg)');
    xval = get(gca,'xtick');
    l = cell(length(xval),1);
    for k = 1:length(l)
        l{k} = sprintf('%.2f',10^xval(k));
    end
    colorbar;
    clim([0, 0.4359]);
    set(gca,'xticklabel',l);
    title(sprintf("Peak Ori: %.0f, Peak SF: %.3f", peak_ori(kern), peak_sf(kern)));
    pause;

end

%%

% figure; hold on;
% sfdata = vertcat(stat_struct.resp_sf);
% sfdata_noise = sfdata(idx_noise, :);
% sfdata = sfdata(idx_noise == 0, :);
% % sfdata = sfdata(vertcat(stat_struct.snr) > 2, :);
% % sfdata = sfdata./max(sfdata, [], 1);
% % sfdata = sfdata .* vertcat(stat_struct.snr);
% plot(sf, sfdata_noise, "-o", "LineWidth", 2, "Color", [0.8, 0.8, 0.8]);
% plot(sf, sfdata, "-o", "LineWidth", 2);
% plot(sf, mean(sfdata, 1), "-ok", "LineWidth", 2);

sf_est = vertcat(stat_struct.sf_est);
figure; hold on;
scatter(1, sf_est(idx_noise));
scatter(2, sf_est(idx_noise == 0));

% figure; hold on;
% scatter(vertcat(stat_struct.sf_est), vertcat(stat_struct.snr));
% 
% figure; hold on;
% scatter(vertcat(stat_struct.sf_est), vertcat(stat_struct.noise));
% 
% figure; hold on;
% histogram(vertcat(stat_struct.tmax))
% 
figure; hold on;
histogram(sf_est(idx_noise == 0));
xline(median(sf_est(idx_noise == 0)), "--r");
title(sprintf("Median sfest = %.3f", median(sf_est(idx_noise == 0))));

%%

% %%%% Analyze full-field drifting gratings
% 
% subj_id = "pv_vip_00";
% recording = "210_000";
% 
% exp_id = char(strcat(subj_id, "_", recording));
% fprintf("Processing %s...\n", exp_id);
% exp_path = char(strcat(exp_id, filesep));
% plane_path = strcat(exp_path, "suite2p", filesep, "plane0", filesep);
% % Generate .align, .segment, and .signals1 files. 
% sbxsuite2sbx(strcat(plane_path, "Fall"), strcat(exp_path, exp_id));
% % Generate .signals file.
% sbxf2spks(strcat(exp_path, exp_id));
% 
% load([strcat(exp_path, exp_id), '.signals'],'-mat');    % load signals




















