function [eye, imgs] = load_eye_data(path, animal, exp, id)

    % Assemble a full or relative path
    if nargin < 4
        fpath = strcat(path, "_", animal, "_", exp);
    else
        fpath = strcat(path, filesep, animal, filesep, animal, "_", exp, "_", id);
    end
    % Check for realtime or suite2p data and load spikes and traces
    d = dir(fpath);
    fnames = {d.name};
    eye_idx = contains(fnames, "_eye.mat");
    if any(eye_idx)
        fname = strcat(fpath, filesep, fnames{eye_idx});
        variable_info = who("-file", fname);
        if ~ismember("eye", variable_info)
            sbxeyemotion(fpath, fname, [2, 50]);
        end
        load(fname, "-mat", "eye", "data");
        imgs = squeeze(data);
    end
end