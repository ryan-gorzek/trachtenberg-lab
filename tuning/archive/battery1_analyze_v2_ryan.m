
% real time analysis...

d = dir('*.sbx'); % get the file name
fname = d(1).name(1:end-4);

%%%% Scanbox metadata
load([fname '.mat']); % contains the syncs and it is loaded into an info variable in the global space

frame_on  = info.evt.stim_on.frame;   % stim onsets
frame_off = info.evt.stim_off.frame;  % stim offsets

%%%% Real-time data
load([fname '_realtime.mat']);        % loads real time signals (raw fluorescence (rtdata frames x ncell matrix) and deconvolved spikes (stdata) )

% Median-filtered traces
rt = rtdata;
% rt = bsxfun(@minus,mean(rtdata),rtdata);    % mean subtract
rt(:,isnan(mean(rtdata))) = [];             % remove bad ROIs
rt = medfilt1(rt, 101) - rt;
rt = zscore(rt);

% Spikes
st = stdata;
st(:,isnan(mean(stdata))) = [];             % remove bad ROIs

ncells = size(rt, 2);

%%%% Stimulus data
load([fname '_p.mat']);              % parameters of stimuli (param) and order of presentation (perms)
                                     % perms is a cell array where each
                                     % entry is the permutation of stimuli
                                     % for that block ...

% this is reconstructing the stim matrix because I didn't save the matrices
% in the _p file.
p1 = unique(param(:,1)); %#ok<*USENS> 
p2 = unique(param(:,2));
nparam1 = length(p1);
nparam2 = length(p2);
[pp1,pp2] = meshgrid(p1,p2);

% Get the stim order
idx = horzcat(perms{:});

%%%% Main code

%%%% Plotting
f = figure; t = tiledlayout(4, 2);
set(f, "Color", "w");

% Plot the responses for traces
nexttile; hold on; plot_stim_resp(rt, frame_on, [-30, 30]); title("Temporal Response, Traces");
% Plot the responses for spikes
nexttile; hold on; plot_stim_resp(st, frame_on, [-30, 30]); title("Temporal Response, Spikes");

% Get the stimulus response matrix for traces
[navg_mean_traces, navg_norm_traces] = tuning(rt, frame_on, frame_off, param, idx, [3, 13], false);
% Get the stimulus response matrix for spikes
[navg_mean_spikes, navg_norm_spikes] = tuning(st, frame_on, frame_off, param, idx);

% average across all orientations to generate a mean spatial frequency tuning curve
sf_mean_traces = squeeze(mean(navg_mean_traces, 2));
sf_norm_traces = mean(navg_norm_traces, 2);
sf_mean_spikes = squeeze(mean(navg_mean_spikes, 2));
sf_norm_spikes = mean(navg_norm_spikes, 2);

% % plot all cells with high sf preference
% idx_high = navg_mean_traces()
% sf_high = navg_mean_traces(:, :, )

% for c = 1:10:ncells
%     figure; imagesc(navg_mean_traces(:, :, c)); colorbar;
%     xticks(1:numel(p1)); xticklabels(round(p1, 3));
%     yticks(1:numel(p2)); yticklabels(round(p2, 3));
%     xlabel("Ori"); ylabel('SF');
% end

%%%% Spatial frequency for all cells, traces
% figure; % tiledlayout(1, 2);
nexttile;
plot(log10(p2), sf_mean_traces, '-o', "LineWidth", 2)   % plot it...
xticks(log10(p2)); xticklabels(round(p2, 3));
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
set(gca, "FontSize", 9, ...
         "XColor", "k", "YColor", "k", ...
         "TickDir", "out", "TickLength", [0.01, 0.01], ...
         "box", "off", ...
         "LineWidth", 1);

%%%% Spatial frequency for all cells, spikes
% figure; % tiledlayout(1, 2);
nexttile;
plot(log10(p2), sf_mean_spikes, '-o', "LineWidth", 2)   % plot it...
xticks(log10(p2)); xticklabels(round(p2, 3));
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
set(gca, "FontSize", 9, ...
         "XColor", "k", "YColor", "k", ...
         "TickDir", "out", "TickLength", [0.01, 0.01], ...
         "box", "off", ...
         "LineWidth", 1);

%%%% Spatial frequency for mean response of all cells, traces
nexttile; hold on;
plot(log10(p2), mean(sf_mean_traces, 2), '-ok', "LineWidth", 2)   % plot it...
xticks(log10(p2)); xticklabels(round(p2, 3));
ylabel('Population Response (Mean)');
[~, x, pv, opt_sf] = fitcurve(log10(p2(1:end-1)), mean(sf_mean_traces(1:end-1, :), 2));
plot(x, pv, "-ob", "LineWidth", 2); xline(log10(opt_sf), "--r", "LineWidth", 2);
title('Population Mean, Traces', ['Optimal sfreq = ' num2str(opt_sf,2)]);  % choose that as the center (optimal) sf
set(gca, "FontSize", 9, ...
         "XColor", "k", "YColor", "k", ...
         "TickDir", "out", "TickLength", [0.01, 0.01], ...
         "box", "off", ...
         "LineWidth", 1);

%%%% Spatial frequency for mean response of all cells, traces
nexttile; hold on;
plot(log10(p2), mean(sf_mean_spikes, 2), '-ok', "LineWidth", 2)   % plot it...
xticks(log10(p2)); xticklabels(round(p2, 3));
ylabel('Population Response (Mean)');
[~, x, pv, opt_sf] = fitcurve(log10(p2(1:end-1)), mean(sf_mean_spikes(1:end-1, :), 2));
plot(x, pv, "-ob", "LineWidth", 2); xline(log10(opt_sf), "--r", "LineWidth", 2);
title('Population Mean, Spikes', ['Optimal sfreq = ' num2str(opt_sf,2)]);  % choose that as the center (optimal) sf
set(gca, "FontSize", 9, ...
         "XColor", "k", "YColor", "k", ...
         "TickDir", "out", "TickLength", [0.01, 0.01], ...
         "box", "off", ...
         "LineWidth", 1);

%%%% Spatial frequency for norm of population response, traces
nexttile; hold on;
plot(log10(p2), sf_norm_traces, '-ok', "LineWidth", 2)   % plot it...
xticks(log10(p2)); xticklabels(round(p2, 3));
xlabel('Log10(spatial freq)');
ylabel('Population Response (Norm)');
[~, x, pv, opt_sf] = fitcurve(log10(p2(1:end-1)), sf_norm_traces(1:end-1, :));
plot(x, pv, "-ob", "LineWidth", 2); xline(log10(opt_sf), "--r", "LineWidth", 2);
title('Population Norm, Traces', ['Optimal sfreq = ' num2str(opt_sf,2)]);  % choose that as the center (optimal) sf
set(gca, "FontSize", 9, ...
         "XColor", "k", "YColor", "k", ...
         "TickDir", "out", "TickLength", [0.01, 0.01], ...
         "box", "off", ...
         "LineWidth", 1);

%%%% Spatial frequency for norm of population response, spikes
nexttile; hold on;
plot(log10(p2), sf_norm_spikes, '-ok', "LineWidth", 2)   % plot it...
xticks(log10(p2)); xticklabels(round(p2, 3));
xlabel('Log10(spatial freq)');
ylabel('Population Response (Norm)');
[~, x, pv, opt_sf] = fitcurve(log10(p2(1:end-1)), sf_norm_spikes(1:end-1, :));
plot(x, pv, "-ob", "LineWidth", 2); xline(log10(opt_sf), "--r", "LineWidth", 2);
title('Population Norm, Spikes', ['Optimal sfreq = ' num2str(opt_sf,2)]);  % choose that as the center (optimal) sf
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
% figure; hold on;
xline(-15, "--k", "LineWidth", 2);
xline(0, "--k", "LineWidth", 2);
xline(15, "--k", "LineWidth", 2);
plot(window(1) : window(2), mean(resp, 3), "LineWidth", 2);
plot(window(1) : window(2), mean(resp, [1, 3]), "LineWidth", 4, "Color", "k");
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
