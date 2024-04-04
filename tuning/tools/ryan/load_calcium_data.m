function [traces, spikes, ops, stat] = load_calcium_data(path, animal, exp, id) % 

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
    % Check for realtime or suite2p data and load spikes and traces
    d = dir(fpath);
    fnames = {d.name};
    suite2p_idx = contains(fnames, "_suite.signals");
    realtime_idx = contains(fnames, "_realtime.mat");
    if any(suite2p_idx)
        load(strcat(fpath, filesep, fnames{suite2p_idx}), "-mat", "sig", "spks");
        traces = sig;
        spikes = spks;
        % load size of imaging field and get spatial location of each cell
        load(strcat(concat_fpath, filesep, "suite2p", filesep, "plane0", filesep, "Fall.mat"), "-mat", "ops", "iscell", "stat");
        stat = stat(iscell(:, 1) > 0);
    elseif any(realtime_idx)
        load(strcat(fpath, filesep, fnames{realtime_idx}), "-mat", "rtdata", "stdata");
        traces = rtdata;
        spikes = stdata;
    end
end