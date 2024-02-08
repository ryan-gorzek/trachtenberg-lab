
%%%% add tools to path
addpath(genpath(fullfile("E:/Imaging/code/tools/")));

subj = "pv_vip_03";
sessions = ["410", "020", "030", "040", "050", "060"];
runs = ["000", "000", "000", "000", "000", "000"];
stims = ["trinoise", "randorisf", "battery1", "battery2", "battery3", "battery4"];
sid = 6;

% load calcium data
[rtdata, stdata, xy, grid] = load_calcium_data(subj, sessions(sid), runs(sid));

% load stimulus data
[frame_on, frame_off, params] = load_stimulus_data(subj, sessions(sid), runs(sid));

% load quadrature data
quad = load_quad_data(subj, sessions(sid), runs(sid));

% load eye data
eye = load_eye_data(subj, sessions(sid), runs(sid));

% create experiment structure
[spikes, dFF] = expstruct(stims(sid), rtdata, stdata, xy, grid, frame_on, frame_off, params, eye, quad);

if stims(sid) == "randorisf", SNR = vertcat(spikes.stats.SNR) > spikes.calc.SNR_thr; end

%% plot center/cross/iso/blank PSTH for context modulation ranked by mean

rank = 2;

stat = vertcat(dFF.stats.cross_mean);
[vals, idx] = sort(stat(:, 2), "descend");
% idx = idx(SNR);
ids = idx(rank);
activity = dFF;

sizes = activity.stimdata.param(activity.stimdata.idx, 5);

data = [];
for st = ["Blank", "Center", "Surround", "Iso", "Cross"] %%%% add surround
    data.(st) = [];
    for id = ids'
        data.(st) = vertcat(data.(st), activity.psths(id).mat(activity.stimdata.context.(lower(st)) & sizes == 10, :));
    end
end

f = figure; tiledlayout(2, 5, "TileSpacing", "tight");
% plot PSTHs
for st = ["Blank", "Center", "Surround", "Iso", "Cross"]
    nexttile; hold on;
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
% scatter responses
for st = {["Blank", "Cross"], ["Center", "Surround"], ["Center", "Iso"], ["Center", "Cross"], ["Cross", "Iso"]}
    nexttile; hold on;
    xline(0, "--k", "LineWidth", 1);
    yline(0, "--k", "LineWidth", 1);
    plot(-1:1, -1:1, "--k", "LineWidth", 1);
    X = vertcat(activity.stats.(lower(st{1}(1)) + "_mean"));
    X = X(SNR, :);
    X(X > 0.5) = 0.5;
    Y = vertcat(activity.stats.(lower(st{1}(2)) + "_mean"));
    Y = Y(SNR, :);
    Y(Y > 0.5) = 0.5;
    C = repmat([0, 0, 0], [nnz(SNR), 1]);
    C(idx == ids, :) = [1, 0, 0];
%     C = SNR;
    scatter(X(:, 2), Y(:, 2), 10, C, "filled", "MarkerFaceAlpha", 0.15, "MarkerEdgeAlpha", 0.15);
    xlabel(st{1}(1) + " [\DeltaF/F]"); ylabel(st{1}(2) + " [\DeltaF/F]");
    % title(st);
    xlim([-0.1, 0.5]); xticks(-0.1:0.1:0.5);
    ylim([-0.1, 0.5]); yticks(-0.1:0.1:0.5);
    axis square;
end

set(f, "Color", "w");

f.Position = [53 252 1388 545];

%% plot mean PSTHs for all cells that pass SNR threshold

rank = 100;

[vals, idx] = sort(vertcat(dFF.stats.cri_mod_mean), "descend");
ids = idx(SNR);
activity = dFF;

pref_con_cent = 0.1;
pref_con_surr = 0.1;

stims = activity.stimdata.param(activity.stimdata.idx, :);

pref_iso = (stims(:, 3) >= pref_con_cent) & ...
           (stims(:, 4) >= pref_con_surr);

pref_cross = (stims(:, 3) >= pref_con_cent) & ...
             (stims(:, 4) >= pref_con_surr);

pref_cb = (stims(:, 3) >= pref_con_cent);

pref_sb = (stims(:, 4) >= pref_con_surr);

data = [];
cells = [];
for st = ["Blank", "Center", "Surround", "Iso", "Cross"]
    cells.(st) = [];
    for id = ids'
        data.(st) = [];
        if st == "Iso"
            data.(st) = vertcat(data.(st), activity.psths(id).mat(activity.stimdata.context.(lower(st)) & pref_iso, :)); % 
        elseif st == "Cross"
            data.(st) = vertcat(data.(st), activity.psths(id).mat(activity.stimdata.context.(lower(st)) & pref_cross, :)); % 
        elseif st == "Center"
            data.(st) = vertcat(data.(st), activity.psths(id).mat(activity.stimdata.context.(lower(st)) & pref_cb, :)); %
        elseif st == "Surround"
            data.(st) = vertcat(data.(st), activity.psths(id).mat(activity.stimdata.context.(lower(st)) & pref_sb, :)); % 
        else
            data.(st) = vertcat(data.(st), activity.psths(id).mat(activity.stimdata.context.(lower(st)), :));
        end
        cells.(st) = vertcat(cells.(st), mean(data.(st), 1));
    end
end

f = figure; tiledlayout(1, 5, "TileSpacing", "tight");
% plot PSTHs
for st = ["Blank", "Center", "Surround", "Iso", "Cross"]
    nexttile; hold on;
    yline(0, "--r", "LineWidth", 1);
    xline(-15, "--k", "LineWidth", 2);
    xline(0, "--k", "LineWidth", 2);
    xline(15, "--k", "LineWidth", 2);
    % plot(activity.psths(1).window(1) : activity.psths(1).window(2), data.(st), "LineWidth", 1, "Color", [0.8, 0.8, 0.8, 0.7]);
    plot(activity.psths(1).window(1) : activity.psths(1).window(2), cells.(st), "LineWidth", 1, "Color", [0, 0, 0, 0.2]);
    title(st);
    axis tight;
    ylim([-0.25, 0.5]);
end

set(f, "Color", "w");

f.Position = [53 252 1388 275];

%% fix a line to the onset of center and surround activity

rank = 100;

[vals, idx] = sort(vertcat(dFF.stats.cri_mod_mean), "descend");
ids = idx(SNR);
activity = dFF;

pref_con_cent = 0.4;
pref_con_surr = 0.4;

stims = activity.stimdata.param(activity.stimdata.idx, :);

pref_iso = (stims(:, 3) >= pref_con_cent) & ...
           (stims(:, 4) >= pref_con_surr);

pref_cross = (stims(:, 3) >= pref_con_cent) & ...
             (stims(:, 4) >= pref_con_surr);

pref_cb = (stims(:, 3) >= pref_con_cent);

pref_sb = (stims(:, 4) >= pref_con_surr);

data = [];
cells = [];
for st = ["Blank", "Center", "Surround", "Iso", "Cross"]
    cells.(st) = [];
    for id = ids'
        data.(st) = [];
        if st == "Iso"
            data.(st) = vertcat(data.(st), activity.psths(id).mat(activity.stimdata.context.(lower(st)) & pref_iso, :)); % 
        elseif st == "Cross"
            data.(st) = vertcat(data.(st), activity.psths(id).mat(activity.stimdata.context.(lower(st)) & pref_cross, :)); % 
        elseif st == "Center"
            data.(st) = vertcat(data.(st), activity.psths(id).mat(activity.stimdata.context.(lower(st)) & pref_cb, :)); %
        elseif st == "Surround"
            data.(st) = vertcat(data.(st), activity.psths(id).mat(activity.stimdata.context.(lower(st)) & pref_sb, :)); % 
        else
            data.(st) = vertcat(data.(st), activity.psths(id).mat(activity.stimdata.context.(lower(st)), :));
        end
        cells.(st) = vertcat(cells.(st), mean(data.(st), 1));
    end
end

wdw = activity.psths(1).window(1) : activity.psths(1).window(2);
slopes = zeros(size(cells.Blank, 1), 2);
for st = ["Center", "Surround"]
    for c = 1:size(cells.(st), 1)
        x = 1:5;
        y = cells.(st)(c, ismember(wdw, x));
        P = polyfit(x, y, 1);
        slopes(c, ["Center", "Surround"] == st) = P(1);
    end
end

figure;
nexttile; hold on;
xline(0, "--k", "LineWidth", 1);
yline(0, "--k", "LineWidth", 1);
plot(0:1, 0:1, "--k", "LineWidth", 1);
X = slopes(:, 1);
Y = slopes(:, 2);
% cmap = jet(9); C = cmap(GRID(SNR), :);
C = GRID(SNR); % repmat([0, 0, 0], [size(Y, 1), 1]);
scatter(X, Y, 20, C, "filled", "MarkerFaceAlpha", 1, "MarkerEdgeAlpha", 1); colorbar;
xline(median(X), "-", "LineWidth", 1, "Color", [1, 0, 0, 0.2]);
yline(median(Y), "-", "LineWidth", 1, "Color", [1, 0, 0, 0.2]);
% if scon == 0.4, xlabel([st{1}(1) + " [\DeltaF/F]", "\bf{" + sprintf("Center Contrast = %i%%", 200*ccon) + "}"]); end
% if ccon == 0.1, ylabel(["{\bf" + sprintf("Surround Contrast = %i%%", 200*scon) + "}", st{1}(2) + " [\DeltaF/F]"]); end
xlim([-0.01, 0.1]); xticks(-0.01:0.05:0.1); xlabel("Center Slope");
ylim([-0.01, 0.1]); yticks(-0.01:0.05:0.1); ylabel("Surround Slope");

set(gcf, "Color", "w");

% f = figure; tiledlayout(1, 5, "TileSpacing", "tight");
% % plot PSTHs
% for st = ["Blank", "Center", "Surround", "Iso", "Cross"]
%     nexttile; hold on;
%     yline(0, "--r", "LineWidth", 1);
%     xline(-15, "--k", "LineWidth", 2);
%     xline(0, "--k", "LineWidth", 2);
%     xline(15, "--k", "LineWidth", 2);
%     % plot(activity.psths(1).window(1) : activity.psths(1).window(2), data.(st), "LineWidth", 1, "Color", [0.8, 0.8, 0.8, 0.7]);
%     plot(activity.psths(1).window(1) : activity.psths(1).window(2), cells.(st), "LineWidth", 1, "Color", [0, 0, 0, 0.2]);
%     title(st);
%     axis tight;
%     ylim([-0.25, 0.5]);
% end
% 
% set(f, "Color", "w");
% 
% f.Position = [53 252 1388 275];

%% plot iso/cross responses at different contrasts

activity = dFF;
ncells = activity.ncells;
wdw = activity.psths(1).window(1) : activity.psths(1).window(2);

stims = activity.stimdata.param(activity.stimdata.idx, :);
f = figure; t = tiledlayout(3, 3, "TileSpacing", "tight");

all_data = {cell(1, 3), cell(1, 3), cell(1, 3)};
for scon = [0.1, 0.2, 0.4]
    for ccon = [0.1, 0.2, 0.4]

        pref.iso = (stims(:, 3) == ccon) & ...
                   (stims(:, 4) == scon);
        
        pref.cross = (stims(:, 3) == ccon) & ...
                     (stims(:, 4) == scon);
        
        pref.center = (stims(:, 3) == ccon);
        
        pref.surround = (stims(:, 4) == scon);
        
        pref.blank = (stims(:, 3) == 0) & ...
                     (stims(:, 4) == 0);
        
        cats = ["Blank", "Center", "Surround", "Iso", "Cross"];
        for st = cats
            data.(st) = zeros(nnz(SNR), 1);
            for id = find(SNR)'
                    celldata = activity.psths(id).mat(activity.stimdata.context.(lower(st)) & pref.(lower(st)), :);
                    data.(st)(find(SNR) == id, 1) = mean(celldata(:, ismember(wdw, 5:20)) - mean(celldata(:, ismember(wdw, -2:3)), 2), "all");
            end
        end

        % scatter responses
        %%%% {["Blank", "Cross"], ["Center", "Surround"], ["Center", "Iso"], ["Center", "Cross"], ["Cross", "Iso"]}
        for st = {["Cross", "Iso"]}
            nexttile; hold on;
            xline(0, "--k", "LineWidth", 1);
            yline(0, "--k", "LineWidth", 1);
            plot(0:1, 0:1, "--k", "LineWidth", 1);
            X = data.(st{1}(1));
            X(X > 0.5) = 0.5;
%             X = X(SNR);
            Y = data.(st{1}(2));
            Y(Y > 0.5) = 0.5;
%             Y = Y(SNR);
            C = repmat([0, 0, 0], [nnz(SNR), 1]);
%             C = C(SNR, :);
            scatter(X, Y, 20, C, "filled", "MarkerFaceAlpha", 0.25, "MarkerEdgeAlpha", 0.25);
            % scatter(median(X), median(Y), "+r");
            xline(median(X), "-", "LineWidth", 1, "Color", [1, 0, 0, 0.2]);
            yline(median(Y), "-", "LineWidth", 1, "Color", [1, 0, 0, 0.2]);
            if scon == 0.4, xlabel([st{1}(1) + " [\DeltaF/F]", "\bf{" + sprintf("Center Contrast = %i%%", 200*ccon) + "}"]); end
            if ccon == 0.1, ylabel(["{\bf" + sprintf("Surround Contrast = %i%%", 200*scon) + "}", st{1}(2) + " [\DeltaF/F]"]); end
            xlim([-0.1, 0.5]); xticks(-0.1:0.1:0.5);
            ylim([-0.1, 0.5]); yticks(-0.1:0.1:0.5);
        end

        all_data{[0.1, 0.2, 0.4] == scon}{[0.1, 0.2, 0.4] == ccon} = data;

    end
end

set(f, "Color", "w");

f.Position = [675.2857 222.7143 974.8571 828.0000];

%% plot center/surround responses at different contrasts

activity = dFF;
ncells = activity.ncells;
wdw = activity.psths(1).window(1) : activity.psths(1).window(2);

stims = activity.stimdata.param(activity.stimdata.idx, :);
f = figure; t = tiledlayout(3, 3, "TileSpacing", "tight");

all_data = {cell(1, 3), cell(1, 3), cell(1, 3)};
for scon = [0.1, 0.2, 0.4]
    for ccon = [0.1, 0.2, 0.4]

        pref.iso = (stims(:, 3) == ccon) & ...
                   (stims(:, 4) == scon);
        
        pref.cross = (stims(:, 3) == ccon) & ...
                     (stims(:, 4) == scon);
        
        pref.center = (stims(:, 3) == ccon);
        
        pref.surround = (stims(:, 4) == scon);
        
        pref.blank = (stims(:, 3) == 0) & ...
                     (stims(:, 4) == 0);
        
        cats = ["Blank", "Center", "Surround", "Iso", "Cross"];
        for st = cats
            data.(st) = zeros(ncells, 1);
            for id = 1:dFF.ncells
                    celldata = activity.psths(id).mat(activity.stimdata.context.(lower(st)) & pref.(lower(st)), :);
                    data.(st)(id, 1) = mean(celldata(:, ismember(wdw, 5:20)) - mean(celldata(:, ismember(wdw, -2:3)), 2), "all");
            end
        end

        % scatter responses
        %%%% {["Blank", "Cross"], ["Center", "Surround"], ["Center", "Iso"], ["Center", "Cross"], ["Cross", "Iso"]}
        for st = {["Center", "Surround"]}
            nexttile; hold on;
            xline(0, "--k", "LineWidth", 1);
            yline(0, "--k", "LineWidth", 1);
            plot(0:1, 0:1, "--k", "LineWidth", 1);
            X = data.(st{1}(1));
            X(X > 0.5) = 0.5;
            X = X(SNR);
            Y = data.(st{1}(2));
            Y(Y > 0.5) = 0.5;
            Y = Y(SNR);
            C = repmat([0, 0, 0], [numel(idx), 1]);
            C = C(SNR, :);
            scatter(X, Y, 20, C, "filled", "MarkerFaceAlpha", 0.25, "MarkerEdgeAlpha", 0.25);
            if scon == 0.4, xlabel([st{1}(1) + " [\DeltaF/F]", "\bf{" + sprintf("Center Contrast = %i%%", 200*ccon) + "}"]); end
            if ccon == 0.1, ylabel(["{\bf" + sprintf("Surround Contrast = %i%%", 200*scon) + "}", st{1}(2) + " [\DeltaF/F]"]); end
            xlim([-0.1, 0.5]); xticks(-0.1:0.1:0.5);
            ylim([-0.1, 0.5]); yticks(-0.1:0.1:0.5);
        end

        all_data{[0.1, 0.2, 0.4] == scon}{[0.1, 0.2, 0.4] == ccon} = data;

    end
end

set(f, "Color", "w");

f.Position = [675.2857 222.7143 974.8571 828.0000];

%%

cons = [0.1, 0.2, 0.4];
figure; tiledlayout(1, 2, "TileSpacing", "tight");
nexttile; hold on;
for scon = 1:3
crossm = []; isom = []; mod_idxm = [];
for ccon = 1:3
    center = all_data{scon}{ccon}.Center;
    cross = all_data{scon}{ccon}.Cross;
    iso = all_data{scon}{ccon}.Iso;
    mod_idx = (cross - iso) ./ (cross + iso);
    crossm = vertcat(crossm, median(cross));
    isom = vertcat(isom, median(iso));
    mod_idxm = vertcat(mod_idxm, median(mod_idx));
end
plot(1:3, mod_idxm, "-", "Color", [1, 1, 1] - cons(scon)*2, "LineWidth", 2);
scatter(1:3, mod_idxm, "o", "filled", "MarkerFaceColor",  [1, 1, 1] - cons(scon)*2);
ylim([-0.1, 1]);
% yyaxis right;
% plot(1:3, crossm, "-b");
% plot(1:3, isom, "-r");
if scon == 1, ylabel("Modulation Index [C - I / C + I]"); end
xticks(1:3); xticklabels(["20%", "40%", "80%"]); xlabel("Center Contrast");
hleg = legend("20%", "", "40%", "", "80%", "", "Box", "off", "Location", "northwest");
title(hleg, "Surround Contrast");
axis square;
end

nexttile; hold on;
for ccon = 1:3
crossm = []; isom = []; mod_idxm = [];
for scon = 1:3
    cross = all_data{scon}{ccon}.Cross;
    iso = all_data{scon}{ccon}.Iso;
    mod_idx = (cross - iso) ./ (cross + iso);
    crossm = vertcat(crossm, median(cross));
    isom = vertcat(isom, median(iso));
    mod_idxm = vertcat(mod_idxm, median(mod_idx));
end
% plot(1:3, crossm, "-b");
% plot(1:3, isom, "-r");
plot(1:3, mod_idxm, "-", "Color", [1, 1, 1] - cons(ccon)*2, "LineWidth", 2);
scatter(1:3, mod_idxm, "o", "filled", "MarkerFaceColor",  [1, 1, 1] - cons(ccon)*2);
ylim([-0.1, 1]);
xticks(1:3); xticklabels(["20%", "40%", "80%"]); xlabel("Surround Contrast");
hleg = legend("20%", "", "40%", "", "80%", "", "Box", "off", "Location", "northwest");
title(hleg, "Center Contrast");
axis square;
end

set(gcf, "Color", "w");

%%

cons = [0.1, 0.2, 0.4];
figure; tiledlayout(1, 2, "TileSpacing", "tight");
nexttile; hold on;
for scon = 1:3
crossm = []; isom = []; mod_idxm = [];
for ccon = 1:3
    center = all_data{scon}{ccon}.Center;
    cross = all_data{scon}{ccon}.Cross;
    iso = all_data{scon}{ccon}.Iso;
    mod_idx = (cross - iso) ./ (cross + iso);
    crossm = vertcat(crossm, median(cross));
    isom = vertcat(isom, median(iso));
    mod_idxm = vertcat(mod_idxm, median(mod_idx));
end
plot(1:3, isom, "-", "Color", [1.2, 0, 0] - [cons(scon)*2, 0, 0], "LineWidth", 2);
scatter(1:3, isom, "o", "filled", "MarkerFaceColor",  [1.2, 0, 0] - [cons(scon)*2, 0, 0]);
ylim([0, 0.2]);
if scon == 1, ylabel("Iso Response [\DeltaF/F]"); end
xticks(1:3); xticklabels(["20%", "40%", "80%"]); xlabel("Center Contrast");
hleg = legend("20%", "", "40%", "", "80%", "", "Box", "off", "Location", "northwest");
title(hleg, "Surround Contrast");
axis square;
end

nexttile; hold on;
for ccon = 1:3
crossm = []; isom = []; mod_idxm = [];
for scon = 1:3
    cross = all_data{scon}{ccon}.Cross;
    iso = all_data{scon}{ccon}.Iso;
    mod_idx = (cross - iso) ./ (cross + iso);
    crossm = vertcat(crossm, median(cross));
    isom = vertcat(isom, median(iso));
    mod_idxm = vertcat(mod_idxm, median(mod_idx));
end
plot(1:3, isom, "-", "Color", [1.2, 0, 0] - [cons(ccon)*2, 0, 0], "LineWidth", 2);
scatter(1:3, isom, "o", "filled", "MarkerFaceColor",  [1.2, 0, 0] - [cons(ccon)*2, 0, 0]);
ylim([0, 0.2]);
xticks(1:3); xticklabels(["20%", "40%", "80%"]); xlabel("Surround Contrast");
hleg = legend("20%", "", "40%", "", "80%", "", "Box", "off", "Location", "northwest");
title(hleg, "Center Contrast");
axis square;
end

set(gcf, "Color", "w");

%%

cons = [0.1, 0.2, 0.4];
figure; tiledlayout(1, 2, "TileSpacing", "tight");
nexttile; hold on;
for scon = 1:3
crossm = []; isom = []; mod_idxm = [];
for ccon = 1:3
    center = all_data{scon}{ccon}.Center;
    cross = all_data{scon}{ccon}.Cross;
    iso = all_data{scon}{ccon}.Iso;
    mod_idx = (cross - iso) ./ (cross + iso);
    crossm = vertcat(crossm, median(cross));
    isom = vertcat(isom, median(iso));
    mod_idxm = vertcat(mod_idxm, median(mod_idx));
end
plot(1:3, crossm, "-", "Color", [0, 0, 1.2, 0] - [0, 0, cons(scon)*2, -cons(scon)*2], "LineWidth", 2);
scatter(1:3, crossm, "o", "filled", "MarkerFaceColor",  [0, 0, 1.2] - [0, 0, cons(scon)*2], "MarkerFaceAlpha", cons(scon)*2);
ylim([0, 0.2]);
if scon == 1, ylabel("Cross Response [\DeltaF/F]"); end
xticks(1:3); xticklabels(["20%", "40%", "80%"]); xlabel("Center Contrast");
hleg = legend("20%", "", "40%", "", "80%", "", "Box", "off", "Location", "southeast");
title(hleg, "Surround Contrast");
axis square;
end

nexttile; hold on;
for ccon = 1:3
crossm = []; isom = []; mod_idxm = [];
for scon = 1:3
    cross = all_data{scon}{ccon}.Cross;
    iso = all_data{scon}{ccon}.Iso;
    mod_idx = (cross - iso) ./ (cross + iso);
    crossm = vertcat(crossm, median(cross));
    isom = vertcat(isom, median(iso));
    mod_idxm = vertcat(mod_idxm, median(mod_idx));
end
plot(1:3, crossm, "-", "Color", [0, 0, 1.2, 0] - [0, 0, cons(ccon)*2, -cons(ccon)*2], "LineWidth", 2);
scatter(1:3, crossm, "o", "filled", "MarkerFaceColor",  [0, 0, 1.2] - [0, 0, cons(ccon)*2], "MarkerFaceAlpha", cons(ccon)*2);
ylim([0, 0.2]);
xticks(1:3); xticklabels(["20%", "40%", "80%"]); xlabel("Surround Contrast");
hleg = legend("20%", "", "40%", "", "80%", "", "Box", "off", "Location", "southeast");
title(hleg, "Center Contrast");
axis square;
end

set(gcf, "Color", "w");

%% plot cross/iso responses at different sizes

activity = dFF;
ncells = activity.ncells;
wdw = activity.psths(1).window(1) : activity.psths(1).window(2);

stims = activity.stimdata.param(activity.stimdata.idx, :);
f = figure; t = tiledlayout(1, 7, "TileSpacing", "tight");

all_data = cell(1, 7);
for sz = [5, 10, 15, 20, 25, 30, 40]

    pref = stims(:, 5) == sz;
    
    cats = ["Blank", "Center", "Surround", "Iso", "Cross"];
    for st = cats
        data.(st) = zeros(ncells, 1);
        for id = 1:dFF.ncells
                celldata = activity.psths(id).mat(activity.stimdata.context.(lower(st)) & pref, :);
                data.(st)(id, 1) = mean(celldata(:, ismember(wdw, 5:20)) - mean(celldata(:, ismember(wdw, -2:3)), 2), "all");
        end
    end

    % scatter responses
    %%%% {["Blank", "Cross"], ["Center", "Surround"], ["Center", "Iso"], ["Center", "Cross"], ["Cross", "Iso"]}
    for st = {["Cross", "Iso"]}
        nexttile; hold on;
        xline(0, "--k", "LineWidth", 1);
        yline(0, "--k", "LineWidth", 1);
        plot(0:1, 0:1, "--k", "LineWidth", 1);
        X = data.(st{1}(1));
        X(X > 0.5) = 0.5;
        X = X(SNR);
        Y = data.(st{1}(2));
        Y(Y > 0.5) = 0.5;
        Y = Y(SNR);
        C = repmat([0, 0, 0], [numel(idx), 1]);
        C = C(SNR, :);
        scatter(X, Y, 20, C, "filled", "MarkerFaceAlpha", 0.25, "MarkerEdgeAlpha", 0.25);
        xline(median(X), "-r"); yline(median(Y), "-r");
        xlabel([st{1}(1) + " [\DeltaF/F]", "\bf{" + sprintf("Center Size = %i%%", sz) + "}"]);
        if sz == 5, ylabel(st{1}(2) + " [\DeltaF/F]"); end
        xlim([-0.1, 0.5]); xticks(-0.1:0.1:0.5);
        ylim([-0.1, 0.5]); yticks(-0.1:0.1:0.5);
    end

    all_data{[5, 10, 15, 20, 25, 30, 40] == sz} = data;

end

set(f, "Color", "w");

f.Position = [231.2857 577 1.7126e+03 353.1429];

%% plot center/cross responses at different sizes

activity = dFF;
ncells = activity.ncells;
wdw = activity.psths(1).window(1) : activity.psths(1).window(2);

stims = activity.stimdata.param(activity.stimdata.idx, :);
f = figure; t = tiledlayout(1, 7, "TileSpacing", "tight");

all_data = cell(1, 7);
for sz = [5, 10, 15, 20, 25, 30, 40]

    pref = stims(:, 5) == sz;
    
    cats = ["Blank", "Center", "Surround", "Iso", "Cross"];
    for st = cats
        data.(st) = zeros(ncells, 1);
        for id = 1:dFF.ncells
                celldata = activity.psths(id).mat(activity.stimdata.context.(lower(st)) & pref, :);
                data.(st)(id, 1) = mean(celldata(:, ismember(wdw, 5:20)) - mean(celldata(:, ismember(wdw, -2:3)), 2), "all");
        end
    end

    % scatter responses
    %%%% {["Blank", "Cross"], ["Center", "Surround"], ["Center", "Iso"], ["Center", "Cross"], ["Cross", "Iso"]}
    for st = {["Center", "Cross"]}
        nexttile; hold on;
        xline(0, "--k", "LineWidth", 1);
        yline(0, "--k", "LineWidth", 1);
        plot(0:1, 0:1, "--k", "LineWidth", 1);
        X = data.(st{1}(1));
        X(X > 0.5) = 0.5;
        X = X(SNR);
        Y = data.(st{1}(2));
        Y(Y > 0.5) = 0.5;
        Y = Y(SNR);
        C = repmat([0, 0, 0], [numel(idx), 1]);
        C = C(SNR, :);
        scatter(X, Y, 20, C, "filled", "MarkerFaceAlpha", 0.25, "MarkerEdgeAlpha", 0.25);
        xlabel([st{1}(1) + " [\DeltaF/F]", "\bf{" + sprintf("Center Size = %i%%", sz) + "}"]);
        if sz == 5, ylabel(st{1}(2) + " [\DeltaF/F]"); end
        xlim([-0.1, 0.5]); xticks(-0.1:0.1:0.5);
        ylim([-0.1, 0.5]); yticks(-0.1:0.1:0.5);
    end

    all_data{[5, 10, 15, 20, 25, 30, 40] == sz} = data;

end

set(f, "Color", "w");

f.Position = [231.2857 577 1.7126e+03 353.1429];

%% plot center/surround responses at different sizes

activity = dFF;
ncells = activity.ncells;
wdw = activity.psths(1).window(1) : activity.psths(1).window(2);

stims = activity.stimdata.param(activity.stimdata.idx, :);
f = figure; t = tiledlayout(1, 7, "TileSpacing", "tight");

all_data = cell(1, 7);
for sz = [5, 10, 15, 20, 25, 30, 40]

    pref = stims(:, 5) == sz;
    
    cats = ["Blank", "Center", "Surround", "Iso", "Cross"];
    for st = cats
        data.(st) = zeros(ncells, 1);
        for id = 1:dFF.ncells
                celldata = activity.psths(id).mat(activity.stimdata.context.(lower(st)) & pref, :);
                data.(st)(id, 1) = mean(celldata(:, ismember(wdw, 5:20)) - mean(celldata(:, ismember(wdw, -2:3)), 2), "all");
        end
    end

    % scatter responses
    %%%% {["Blank", "Cross"], ["Center", "Surround"], ["Center", "Iso"], ["Center", "Cross"], ["Cross", "Iso"]}
    for st = {["Center", "Surround"]}
        nexttile; hold on;
        xline(0, "--k", "LineWidth", 1);
        yline(0, "--k", "LineWidth", 1);
        plot(0:1, 0:1, "--k", "LineWidth", 1);
        X = data.(st{1}(1));
        X(X > 0.5) = 0.5;
        X = X(SNR);
        Y = data.(st{1}(2));
        Y(Y > 0.5) = 0.5;
        Y = Y(SNR);
        C = repmat([0, 0, 0], [numel(idx), 1]);
        C = C(SNR, :);
        scatter(X, Y, 20, C, "filled", "MarkerFaceAlpha", 0.25, "MarkerEdgeAlpha", 0.25);
        xline(median(X), "-r"); yline(median(Y), "-r");
        xlabel([st{1}(1) + " [\DeltaF/F]", "\bf{" + sprintf("Center Size = %i%%", sz) + "}"]);
        if sz == 5, ylabel(st{1}(2) + " [\DeltaF/F]"); end
        xlim([-0.1, 0.5]); xticks(-0.1:0.1:0.5);
        ylim([-0.1, 0.5]); yticks(-0.1:0.1:0.5);
    end

    all_data{[5, 10, 15, 20, 25, 30, 40] == sz} = data;

end

set(f, "Color", "w");

f.Position = [231.2857 577 1.7126e+03 353.1429];

%% plot center/surround responses at different sizes (spikes)

activity = spikes;
ncells = activity.ncells;
wdw = activity.psths(1).window(1) : activity.psths(1).window(2);

stims = activity.stimdata.param(activity.stimdata.idx, :);
f = figure; t = tiledlayout(1, 7, "TileSpacing", "tight");

all_data = cell(1, 7);
for sz = [5, 10, 15, 20, 25, 30, 40]

    pref = stims(:, 5) == sz;
    
    cats = ["Blank", "Center", "Surround", "Iso", "Cross"];
    for st = cats
        data.(st) = zeros(ncells, 1);
        for id = 1:dFF.ncells
                celldata = activity.psths(id).mat(activity.stimdata.context.(lower(st)) & pref, :);
                data.(st)(id, 1) = mean(celldata(:, ismember(wdw, 0:10)), "all");
        end
    end

    % scatter responses
    %%%% {["Blank", "Cross"], ["Center", "Surround"], ["Center", "Iso"], ["Center", "Cross"], ["Cross", "Iso"]}
    for st = {["Center", "Surround"]}
        nexttile; hold on;
        xline(0, "--k", "LineWidth", 1);
        yline(0, "--k", "LineWidth", 1);
        plot(0:1, 0:1, "--k", "LineWidth", 1);
        X = data.(st{1}(1));
        X(X > 0.5) = 0.5;
        X = X(SNR);
        Y = data.(st{1}(2));
        Y(Y > 0.5) = 0.5;
        Y = Y(SNR);
        C = repmat([0, 0, 0], [numel(idx), 1]);
        C = C(SNR, :);
        scatter(X, Y, 20, C, "filled", "MarkerFaceAlpha", 0.25, "MarkerEdgeAlpha", 0.25);
        xlabel([st{1}(1) + " [\DeltaF/F]", "\bf{" + sprintf("Center Size = %i%%", sz) + "}"]);
        if sz == 5, ylabel(st{1}(2) + " [\DeltaF/F]"); end
        xlim([-0.05, 0.1]); xticks(-0.1:0.1:0.5);
        ylim([-0.05, 0.1]); yticks(-0.1:0.1:0.5);
    end

    all_data{[5, 10, 15, 20, 25, 30, 40] == sz} = data;

end

set(f, "Color", "w");

f.Position = [231.2857 577 1.7126e+03 353.1429];

%% boxplot surround distribution versus full-field

activity = dFF;
ncells = activity.ncells;
wdw = activity.psths(1).window(1) : activity.psths(1).window(2);

stims = activity.stimdata.param(activity.stimdata.idx, :);
f = figure; t = tiledlayout(1, 2, "TileSpacing", "tight");

cent_data = [];
surr_data = [];
for sz = [5, 15, 25, 35, 45]

    pref = stims(:, 5) == sz;
    
    cats = ["Blank", "Center", "Surround", "Iso", "Cross"];
    for st = cats
        data.(st) = zeros(ncells, 1);
        for id = 1:dFF.ncells
                celldata = activity.psths(id).mat(activity.stimdata.context.(lower(st)) & pref, :);
                data.(st)(id, 1) = mean(celldata(:, ismember(wdw, 5:20)) - mean(celldata(:, ismember(wdw, -2:3)), 2), "all");
        end
    end

    % gather responses
    if sz == 5, surr_data = horzcat(surr_data, data.Surround); end
    surr_data = horzcat(surr_data, data.Iso);
    cent_data = horzcat(cent_data, data.Iso);
    if sz == 45, cent_data = horzcat(cent_data, data.Center); end

end

nexttile; boxPlot(surr_data);
nexttile; boxPlot(circshift(cent_data, 1, 2));

set(f, "Color", "w");

f.Position = [231.2857 577 1.7126e+03 353.1429];
