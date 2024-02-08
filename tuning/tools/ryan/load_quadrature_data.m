function quad = load_quad_data(path, animal, exp, id)

    % Assemble a full or relative path
    if nargin < 4
        fpath = strcat(path, "_", animal, "_", exp);
    else
        fpath = strcat(path, filesep, animal, "_", exp, "_", id);
    end
    % Check for realtime or suite2p data and load spikes and traces
    d = dir(fpath);
    fnames = {d.name};
    quadrature_idx = contains(fnames, "_quadrature.mat");
    if any(quadrature_idx)
        load(strcat(fpath, filesep, fnames{quadrature_idx}), "-mat", "quad_data");
        quad = abs(diff(quad_data));
    end
end