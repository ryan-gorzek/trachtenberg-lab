
% add tools to path
addpath(genpath(fullfile("E:/Imaging/code/tools/")));

% load calcium data (from realtime for now)
[rtdata, stdata] = load_calcium_data("pv_vip_00", "320", "000");

% load stimulus data (2D for now)
[frame_on, frame_off, params] = load_stimulus_data("pv_vip_00", "320", "000");

%%%% Process calcium data
% Median-filtered traces
rt = rtdata;
rt(:,isnan(mean(rtdata))) = [];             % remove bad ROIs
% delta F / F
nstim = numel(frame_on);
ncells = size(rt, 2);
baseline = zeros(nstim, ncells);
for s = 1:nstim
    baseline(s, :) = mean(rt(frame_on(s) - 2 : frame_on(s) + 3, :), 1);
    % or median of entire recording?
end
rt = (mean(baseline, 1) - rt) ./ mean(baseline, 1);

% Spikes
st = stdata;
st(:,isnan(mean(stdata))) = [];             % remove bad ROIs

% Get the stimulus response matrix for traces
[navg_mean_traces, navg_norm_traces] = tuning(rt, frame_on, frame_off, params.param, params.idx, [3, 3]);
% Get the stimulus response matrix for spikes
[navg_mean_spikes, navg_norm_spikes] = tuning(st, frame_on, frame_off, params.param, params.idx, [0, 0]);

%%%%

params.p2 = 60 ./ params.p2;
p1 = params.p1;
p2 = params.p2;
pp2 = flip(p2);

%% 

[~, peak_idx_spikes] = sort(squeeze(max(navg_mean_traces, [], [1, 2])), "descend");

% for cell_id = peak_idx_spikes(1:50)'

cell_id = peak_idx_spikes(1);

f = figure; set(f, "Color", "w"); tiledlayout(2, 4);
% traces
dFF = rt(:, cell_id);
nexttile([1, 2]); plot(dFF, "LineWidth", 0.5, "Color", "k"); title("\DeltaF / F"); axis tight;
nexttile; hold on; plot_stim_resp(dFF, frame_on, [-45, 45]);
[navg_dFF, ~] = tuning(dFF, frame_on, frame_off, params.param, params.idx, [3, 3], false);
nexttile; imagesc(flip(navg_dFF, 1)); colorbar; clim([-0.05, 0.4]);
xticks(1:numel(params.p1)); xticklabels(round(params.p1, 3)); xlabel(sprintf("Direction (%s)", char(176)));
yticks(1:numel(params.p2)); yticklabels(round(pp2, 3)); ylabel("Temporal Frequency (Hz)");
% spikes
spiking = st(:, cell_id);
nexttile([1, 2]); plot(spiking, "LineWidth", 0.5, "Color", "k"); title("Inferred Spikes"); axis tight;
nexttile; hold on; plot_stim_resp(spiking, frame_on, [-45, 45]);
[navg_spiking, ~] = tuning(spiking, frame_on, frame_off, params.param, params.idx, [0, 0], false);
nexttile; imagesc(flip(navg_spiking, 1), [0, 0.75]); colorbar; clim([0, 0.35]);
xticks(1:numel(params.p1)); xticklabels(round(params.p1, 3)); xlabel(sprintf("Direction (%s)", char(176)));
yticks(1:numel(params.p2)); yticklabels(round(pp2, 3)); ylabel("Temporal Frequency (Hz)");

f.Position = [802.1429 487.2857 1.3577e+03 632];

% pause;
% 
% end

%%

% % rt = bsxfun(@minus, mean(rt), rt);    % mean subtract
% rt(:,isnan(mean(rtdata))) = [];               % remove bad ROIs
% rt = medfilt1(rt, 101) - rt;
% % rt = zscore(rt);
% 
% % Spikes
% st = stdata;
% st(:,isnan(mean(stdata))) = [];             % remove bad ROIs
% 
% ncells = size(rt, 2);



%%%% Main code

%%%% Plotting
f = figure; t = tiledlayout(4, 2);
set(f, "Color", "w");

% Plot the responses for traces
nexttile; hold on; plot_stim_resp(rt, frame_on, [-45, 45]); title("Temporal Response, Traces");
% Plot the responses for spikes
nexttile; hold on; plot_stim_resp(st, frame_on, [-45, 45]); title("Temporal Response, Spikes");

% average across all orientations to generate a mean spatial frequency tuning curve
sf_mean_traces = squeeze(mean(navg_mean_traces, 2));
sf_norm_traces = mean(navg_norm_traces, 2);
sf_mean_spikes = squeeze(mean(navg_mean_spikes, 2));
sf_norm_spikes = mean(navg_norm_spikes, 2);

% % plot all cells with high sf preference
% idx_high = navg_mean_traces()
% sf_high = navg_mean_traces(:, :, )

p2p = flip(p2); %%%%

%%%% Spatial frequency for all cells, traces
% figure; % tiledlayout(1, 2);
nexttile;
plot(log10(p2), sf_mean_traces, '-o', "LineWidth", 2)   % plot it...
xticks(log10(p2p)); xticklabels(round(p2p, 2));
ylabel('Mean Response');
% opt_sf_cells = zeros(ncells, 1);
% for c = 1:ncells
%     p=polyfit(log10(p2), sfcells(:, c), 2);         % fit a parabola (poly deg 2)
%     opt_sf_cells(c) = 10^(-p(2)/(2*p(1)));           % -b/(2a)
% end
% nexttile;
% histogram(log10(opt_sf_cells(opt_sf_cells <= max(p2))));
% xticks(log10(p2)); xticklabels(round(p2, 2));
% ylabel("# of Cells");
axis tight;
set(gca, "FontSize", 9, ...
         "XColor", "k", "YColor", "k", ...
         "TickDir", "out", "TickLength", [0.01, 0.01], ...
         "box", "off", ...
         "LineWidth", 1);

%%%% Spatial frequency for all cells, spikes
% figure; % tiledlayout(1, 2);
nexttile;
plot(log10(p2), sf_mean_spikes, '-o', "LineWidth", 2)   % plot it...
xticks(log10(p2p)); xticklabels(round(p2p, 2));
ylabel('Mean Response');
% opt_sf_cells = zeros(ncells, 1);
% for c = 1:ncells
%     p=polyfit(log10(p2), sf_mean_spikes(:, c), 2);         % fit a parabola (poly deg 2)
%     opt_sf_cells(c) = 10^(-p(2)/(2*p(1)));           % -b/(2a)
% end
% nexttile;
% histogram(log10(opt_sf_cells(opt_sf_cells <= max(p2))));
% xticks(log10(p2)); xticklabels(round(p2, 2));
% ylabel("# of Cells");
axis tight;
set(gca, "FontSize", 9, ...
         "XColor", "k", "YColor", "k", ...
         "TickDir", "out", "TickLength", [0.01, 0.01], ...
         "box", "off", ...
         "LineWidth", 1);

%%%% Spatial frequency for mean response of all cells, traces
nexttile; hold on;
plot(log10(p2), mean(sf_mean_traces, 2), '-ok', "LineWidth", 2)   % plot it...
xticks(log10(p2p)); xticklabels(round(p2p, 2));
ylabel('Population Response (Mean)');
[~, x, pv, opt_sf] = fitcurve(log10(p2), mean(sf_mean_traces, 2));
plot(x, pv, "-ob", "LineWidth", 2); % xline(opt_sf, "--r", "LineWidth", 2);
axis tight;
title('Population Mean, Traces', ['Optimal tfreq = ' num2str(opt_sf,2)]);  % choose that as the center (optimal) sf
set(gca, "FontSize", 9, ...
         "XColor", "k", "YColor", "k", ...
         "TickDir", "out", "TickLength", [0.01, 0.01], ...
         "box", "off", ...
         "LineWidth", 1);

%%%% Spatial frequency for mean response of all cells, traces
nexttile; hold on;
plot(log10(p2), mean(sf_mean_spikes, 2), '-ok', "LineWidth", 2)   % plot it...
xticks(log10(p2p)); xticklabels(round(p2p, 2));
ylabel('Population Response (Mean)');
[~, x, pv, opt_sf] = fitcurve(log10(p2), mean(sf_mean_spikes, 2));
plot(x, pv, "-ob", "LineWidth", 2); % xline(opt_sf, "--r", "LineWidth", 2);
axis tight;
title('Population Mean, Spikes', ['Optimal tfreq = ' num2str(opt_sf,2)]);  % choose that as the center (optimal) sf
set(gca, "FontSize", 9, ...
         "XColor", "k", "YColor", "k", ...
         "TickDir", "out", "TickLength", [0.01, 0.01], ...
         "box", "off", ...
         "LineWidth", 1);

%%%% Spatial frequency for norm of population response, traces
nexttile; hold on;
plot(log10(p2), sf_norm_traces, '-ok', "LineWidth", 2)   % plot it...
xticks(log10(p2p)); xticklabels(round(p2p, 2));
xlabel('Log10(temporal freq)');
ylabel('Population Response (Norm)');
[~, x, pv, opt_sf] = fitcurve(log10(p2), sf_norm_traces);
plot(x, pv, "-ob", "LineWidth", 2); % xline(opt_sf, "--r", "LineWidth", 2);
axis tight;
title('Population Norm, Traces', ['Optimal tfreq = ' num2str(opt_sf,2)]);  % choose that as the center (optimal) sf
set(gca, "FontSize", 9, ...
         "XColor", "k", "YColor", "k", ...
         "TickDir", "out", "TickLength", [0.01, 0.01], ...
         "box", "off", ...
         "LineWidth", 1);

%%%% Spatial frequency for norm of population response, spikes
nexttile; hold on;
plot(log10(p2), sf_norm_spikes, '-ok', "LineWidth", 2)   % plot it...
xticks(log10(p2p)); xticklabels(round(p2p, 2));
xlabel('Log10(temporal freq)');
ylabel('Population Response (Norm)');
[~, x, pv, opt_sf] = fitcurve(log10(p2), sf_norm_spikes);
plot(x, pv, "-ob", "LineWidth", 2); % xline(opt_sf, "--r", "LineWidth", 2);
axis tight;
title('Population Norm, Spikes', ['Optimal tfreq = ' num2str(opt_sf,2)]);  % choose that as the center (optimal) sf
set(gca, "FontSize", 9, ...
         "XColor", "k", "YColor", "k", ...
         "TickDir", "out", "TickLength", [0.01, 0.01], ...
         "box", "off", ...
         "LineWidth", 1);

% % lets look at the temporal response for the spikes
% 
% savg = zeros(1,45);
% for i=1:nstim
%     savg = savg + vecnorm(rt(frame_on(i)-15 : frame_on(i)+29, :),2,2)';
% end

f.Position = [1.1879e+03 41.5714 745.7143 1.1154e+03];

function plot_stim_resp(pop, frameon, window)

nstim = numel(frameon);
ncells = size(pop, 2);
nframes = numel(window(1) : window(2));
resp = zeros(ncells, nframes, nstim);
for n = 1:nstim
    cw = frameon(n) + window(1) : frameon(n) + window(2);
    resp(:, :, n) = pop(cw, :)';
end
if ncells == 1
    resp = permute(resp, [3, 2, 1]);
end
xline(-15, "--k", "LineWidth", 2);
xline(0, "--k", "LineWidth", 2);
xline(30, "--k", "LineWidth", 2);
if ncells == 1
    plot(window(1) : window(2), mean(resp, 3), "LineWidth", 1, "Color", [0.8, 0.8, 0.8, 0.7]);
    plot(window(1) : window(2), std(squeeze(resp)), "LineWidth", 4, "Color", "k");
else
    plot(window(1) : window(2), mean(resp, 3), "LineWidth", 2);
    plot(window(1) : window(2), mean(resp, [1, 3]), "LineWidth", 4, "Color", "k");
end
% set(gcf, "Color", "w");
axis tight;
xlabel("Frames from Stimulus Onset");
set(gca, "FontSize", 9, ...
         "XColor", "k", "YColor", "k", ...
         "TickDir", "out", "TickLength", [0.01, 0.01], ...
         "box", "off", ...
         "LineWidth", 1);
end

function [navg_mean, navg_norm] = tuning(data, frameon, frameoff, param, idx, window, subtract)

if nargin < 6, window = [0, 0]; end
if nargin < 7, subtract = false; end

nstim = numel(frameon);
ncells = size(data, 2);

vmean = zeros(nstim, ncells);
vnorm = zeros(nstim, 1);
for i = 1:nstim                                         % for each stim index
    v = mean(data(frameon(i) + window(1) : frameoff(i) + window(2), :));
    vn = v;
    if subtract == true
        v = v - mean(data(frameon(i) - 4 : frameon(i) - 1, :));
    end
    vmean(i,:) = v;
    vnorm(i) = mean(vecnorm(vn, 2, 2));
end

M = max(idx);
p1 = unique(param(:,1)); %#ok<*USENS> 
p2 = unique(param(:,2));
[pp1,pp2] = meshgrid(p1,p2);
navg_mean = zeros(M, ncells);                           % and the average (across repeats) response to each
navg_norm = zeros(M, 1);
for k = 1:M
    j = idx == k;                                       % find all the times we presented stim k
    navg_mean(k,:) = mean(vmean(j,:));                  % average the response
    navg_norm(k) = mean(vnorm(j));
end

navg_mean = reshape(navg_mean, [size(pp1) ncells]);     % reshape so it has the right matrix form
navg_norm = reshape(navg_norm, size(pp1));

end

function [p, x, pv, opt] = fitcurve(x, y)

p = polyfit(x, y, 2);         % fit a parabola (poly deg 2)
opt = 10^(-p(2)/(2*p(1)));    % -b/(2a)
pv = polyval(p, x);
end
