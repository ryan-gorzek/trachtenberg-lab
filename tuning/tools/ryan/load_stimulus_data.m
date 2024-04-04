function [frameon, frameoff, params, stim_center] = load_stimulus_data(path, animal, exp, id)

    % Assemble a full or relative path
    if nargin < 4
        [animal, exp, id] = deal(path, animal, exp);
        fpath = strcat(animal, "_", exp, "_", id);
        concat_code = split(join(repmat(exp, [1, 3]), ""), "");
        concat_fpath = strcat(animal, "_", join(concat_code(2:3:10), ""), "_", join(concat_code(2:3:10), ""));
    else
        fpath = strcat(path, filesep, animal, filesep, animal, "_", exp, "_", id);
        concat_code = split(join(repmat(exp, [1, 3]), ""), "");
        concat_fpath = strcat(path, filesep, animal, filesep, animal, "_", join(concat_code(2:3:10), ""), "_", join(concat_code(2:3:10), ""));
    end
    % Load timestamps from scanbox metadata
    [~, fname, ~] = fileparts(fpath);
    load(strcat(fpath, filesep, fname, ".mat"), "info");
    % load stimulus center
    stim_center = readmatrix(strcat(concat_fpath, filesep, "suite2p", filesep, "plane0", filesep, "stim_center.csv"));
    % get frameon and frameoff
%     frameon = info.evt.photodiode_on.frame(1:2:end);
%     if any(diff(frameon) < 0)
%         idx = find(diff(frameon) < 0) + 1;
%         frameon(idx : end) = round(frameon(idx : end) + 2^16);
%     end
    frameon = info.evt.stim_on.frame;
    if any(diff(frameon) < 0)
        idx = find(diff(frameon) < 0) + 1;
        frameon(idx : end) = round(frameon(idx : end) + 2^16);
    end
%     frameoff = info.evt.photodiode_on.frame(2:2:end);
    % Check for param or log data and load spikes and traces
    d = dir(fpath);
    fnames = {d.name};
    orisflog_idx = contains(fnames, ".log_31");
    trinoiselog_idx = contains(fnames, ".log_59");
    param_idx = contains(fnames, "_p.mat");
    if any(orisflog_idx) % read orisf data
        log = sbxreadorisflog(char(strcat(fpath, filesep, erase(fnames{orisflog_idx}, ".log_31"))));
        tbl = log{1};
        stims = table2array(tbl(:, ["sper", "ori", "sphase"]));
        combos = unique(stims, "rows");
        oris = 0:10:170;
        sfs = 1./ ((1920 ./ (1.31 .^ (0:11))) ./ 14.3258);
        sphases = linspace(0, 360, 9); sphases = sphases(1:end-1);
        param = horzcat(sfs(combos(:, 1) + 1)', oris(combos(:, 2) + 1)', sphases(combos(:, 3) + 1)');
        params = build_params(param);
        mapping = dictionary(num2cell(combos, 2)', 1:size(combos, 1));
        idx = nan(1, size(stims, 1));
        for st = 1:size(stims, 1)
            idx(st) = mapping({stims(st, :)});
        end
        params.idx = idx;
        params.stimidx = stims + 1;
        % frameon and frameoff for orisf
        frameon = table2array(tbl(:, "sbxframe"));
        frameoff = vertcat(frameon(2:end) - 1, info.evt.stim_off.frame);
    elseif any(trinoiselog_idx)
        trinoisedata_idx = contains(fnames, ".trinoise_kern");
        if ~any(trinoisedata_idx)
            r = sbxtrinoise_psth(char(strcat(fpath, filesep, erase(fnames{trinoiselog_idx}, ".log_59"))), info);
        else
            load(strcat(fpath, filesep, fnames{trinoisedata_idx}), "-mat", "r");
        end
        frameon = []; frameoff = [];
        params = r;
    elseif any(param_idx)
        load(strcat(fpath, filesep, fnames{param_idx}), "-mat", "param", "perms");
        params = build_params(param);
        params.perms = perms;
        params.idx = horzcat(perms{:});
        params.stimidx = get_stimidx(params);
        frameoff = [];
    end

    function params = build_params(param)    
        params.param = param;
        switch size(param, 2)
            case 2
                params.p1 = unique(param(:, 1)); %#ok<*USENS>
                params.p2 = unique(param(:, 2));
                params.np1 = numel(params.p1);
                params.np2 = numel(params.p2);
                [params.pp1, params.pp2] = meshgrid(params.p1, params.p2);
            case 3
                params.p1 = unique(param(:, 1)); %#ok<*USENS>
                params.p2 = unique(param(:, 2));
                params.p3 = unique(param(:, 3));
                params.np1 = numel(params.p1);
                params.np2 = numel(params.p2);
                params.np3 = numel(params.p3);
                [params.pp1, params.pp2, params.pp3] = meshgrid(params.p1, params.p2, params.p3);
            case 4
                params.p1 = unique(param(:, 1)); %#ok<*USENS>
                params.p2 = unique(param(:, 2));
                params.p3 = unique(param(:, 3));
                params.p4 = unique(param(:, 4));
                params.np1 = numel(params.p1);
                params.np2 = numel(params.p2);
                params.np3 = numel(params.p3);
                params.np4 = numel(params.p4);
                [params.pp1, params.pp2, params.pp3, params.pp4] = ndgrid(params.p1, params.p2, params.p3, params.p4);
            case 5
                params.p1 = unique(param(:, 1)); %#ok<*USENS>
                params.p2 = unique(param(:, 2));
                params.p3 = unique(param(:, 3));
                params.p4 = unique(param(:, 4));
                params.p5 = unique(param(:, 5));
                params.np1 = numel(params.p1);
                params.np2 = numel(params.p2);
                params.np3 = numel(params.p3);
                params.np4 = numel(params.p4);
                params.np5 = numel(params.p5);
                [params.pp1, params.pp2, params.pp3, params.pp4, params.pp5] = ndgrid(params.p1, params.p2, params.p3, params.p4, params.p5);
        end
    end

    function stimidx = get_stimidx(params)
        stimvals = params.param(params.idx, :);
        stimidx = nan(size(stimvals));
        for s = 1:size(params.param, 2)
            [~, stimidx(:, s)] = ismember(stimvals(:, s), params.("p" + string(s)));
        end
    end

end