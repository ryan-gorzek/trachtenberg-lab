function [eye, imgs] = load_eye_data(path, animal, exp, id)

    % Assemble a full or relative path
    if nargin < 4
        [animal, exp, id] = deal(path, animal, exp);
        fpath = strcat(animal, "_", exp, "_", id);
    else
        fpath = strcat(path, filesep, animal, filesep, animal, "_", exp, "_", id);
    end
    % get suite2p path for eye center
    concat_code = split(join(repmat(exp, [1, 3]), ""), "");
    concat_fpath = strcat(animal, "_", join(concat_code(2:3:10), ""), "_", join(concat_code(2:3:10), ""));
    % Check for realtime or suite2p data and load spikes and traces
    d = dir(fpath);
    fnames = {d.name};
    eye_idx = contains(fnames, "_eye.mat");
    if any(eye_idx)
        fname = strcat(fpath, filesep, fnames{eye_idx});
        variable_info = who("-file", fname);
        if ~ismember("eye", variable_info)
            sbxeyemotion(fname, [2, 50], concat_fpath);
        end
        load(fname, "-mat", "eye", "data");
        imgs = squeeze(data);
    end
end