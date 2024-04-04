function [spikes, dFF] = expstruct(stim_type, traces, spikes, ops, stat, frameon, frameoff, params, stim_center, pix_per_deg, eye, imgs, eye_center, eye_radius, quad)

%%%% add metadata
expst = [];
expst.ncells = size(traces, 2);
expst.nframes = size(traces, 1);
expst.framerate = 15.62;
if stim_type ~= "trinoise"
    expst.ndims = size(params.param, 2);
    expst.nstims = numel(frameon);
    expst.stimshape = size(params.pp1);
    expst.stimdata = params;
else
    imsize = [512, 796];
    expst.maxproj = zeros([imsize, 3]);
    expst.maxproj(ops.yrange(1) : ops.yrange(2)-1, ops.xrange(1) : ops.xrange(2)-1, 1) = ops.max_proj;
    expst.maxproj = uint8((expst.maxproj ./ max(expst.maxproj, [], "all")) .* 255);
    expst.masks = zeros([imsize, expst.ncells]);
    for c = 1:expst.ncells
        mask = zeros(imsize);
        ind = sub2ind(size(mask), stat{c}.ypix, stat{c}.xpix);
        mask(ind) = 1;
        expst.masks(:, :, c) = mask;
    end
end

switch stim_type
    case "trinoise"
        data.spikes = expst; data.dFF = expst;
        data.spikes.stats = params; data.dFF.stats = params;
        ncells = numel(data.spikes.stats);
        stim_center(2) = 1080 - stim_center(2);
        data.spikes.stim_center = stim_center; data.dFF.stim_center = stim_center;
        %%%% calculate distance from stimulus center
        for activity = ["spikes", "dFF"]
            on_dist = cell(1, ncells); off_dist = cell(1, ncells);
            on_dist_fit = cell(1, ncells); off_dist_fit = cell(1, ncells);
            for c = 1:ncells
                ON = 80 .* [data.(activity).stats(c).on_x + 1, data.(activity).stats(c).on_y + 1];
                OFF = 80 .* [data.(activity).stats(c).off_x + 1, data.(activity).stats(c).off_y + 1];
                ON_fit = 80 .* [data.(activity).stats(c).on_x_fit + 1, data.(activity).stats(c).on_y_fit + 1];
                OFF_fit = 80 .* [data.(activity).stats(c).off_x_fit + 1, data.(activity).stats(c).off_y_fit + 1];
                on_dist{c} = sqrt(sum((ON - stim_center) .^ 2)) / pix_per_deg;
                off_dist{c} = sqrt(sum((OFF - stim_center) .^ 2)) / pix_per_deg;
                on_dist_fit{c} = sqrt(sum((ON_fit - stim_center) .^ 2)) / pix_per_deg;
                off_dist_fit{c} = sqrt(sum((OFF_fit - stim_center) .^ 2)) / pix_per_deg;
                % recalculate SNR
                data.(activity).stats(c).snr_on_all = data.(activity).stats(c).snr_on;
                data.(activity).stats(c).snr_on = data.(activity).stats(c).snr_on(data.(activity).stats(c).tmax_on);
                if isinf(data.(activity).stats(c).snr_on), data.(activity).stats(c).snr_on = NaN; end
                data.(activity).stats(c).snr_off_all = data.(activity).stats(c).snr_off;
                data.(activity).stats(c).snr_off = data.(activity).stats(c).snr_off(data.(activity).stats(c).tmax_off);
                if isinf(data.(activity).stats(c).snr_off), data.(activity).stats(c).snr_off = NaN; end
            end
            [data.(activity).stats.on_dist] = on_dist{:}; [data.(activity).stats.off_dist] = off_dist{:};
            [data.(activity).stats.on_dist_fit] = on_dist_fit{:}; [data.(activity).stats.off_dist_fit] = off_dist_fit{:};
            % calculate SNR thresholds
            data.(activity).calc.snr_on_thr = SNR_thr(vertcat(data.(activity).stats.snr_on), vertcat(data.(activity).stats.tmax_on));
            data.(activity).calc.snr_off_thr = SNR_thr(vertcat(data.(activity).stats.snr_off), vertcat(data.(activity).stats.tmax_off));
            % add eye data
            data.(activity).eye.img = imgs(:, :, end);
            eye_data = 2 .* vertcat(eye.Radius); % use diameter
            data.(activity).eye.data = medfilt1(eye_data(1:expst.nframes), 15, "truncate"); % median filter to remove blinking artifacts
            data.(activity).eye.pos = vertcat(eye.Centroid);
            data.(activity).eye.pos = data.(activity).eye.pos(1:data.(activity).nframes, :);
            eye_pos_mode = mode(data.(activity).eye.pos, 1);
            data.(activity).eye.center = eye_pos_mode;
            eye_pos_logical = (data.(activity).eye.pos(:, 1) > eye_pos_mode(1) - eye_radius) & ...
                              (data.(activity).eye.pos(:, 1) < eye_pos_mode(1) + eye_radius) & ...
                              (data.(activity).eye.pos(:, 2) > eye_pos_mode(2) - eye_radius) & ...
                              (data.(activity).eye.pos(:, 2) < eye_pos_mode(2) + eye_radius);
            data.(activity).eye.logical = eye_pos_logical;
        end
        spikes = data.spikes; dFF = data.dFF;
        return
    case "randorisf"
        expst.stimlen = 16;
        %%%% ignore the spatial phase dimension
        expst.ndims = 2;
        expst.stimshape = [18, 12];
        %%%% assign stimulus names and values
        expst.stimabbrv = "ROSF";
        expst.dimnames = [sprintf("Orientation (%s)", char(176)), "Spatial Frequency (cpd)"];
        expst.dimabbrv = ["Ori", "SF"];
        expst.dimvals = {params.p2', params.p1'};
        expst.dimscales = ["polar", "log"];
        %%%% function handles for estimating peak parameters
        expst.calc.estfunc = {@ori_est, @sf_est};
        %%%% add visualization parameters
        expst.viz.window = [-30, 30];
        expst.viz.labels = {sprintf("Orientation (%s)", char(176)), "Spatial Frequency (cpd)"};
        expst.viz.ticks.cont = {expst.stimdata.p2, log10(expst.stimdata.p1)};
        expst.viz.ticks.disc = {1:expst.stimdata.np2, 1:expst.stimdata.np1};
        expst.viz.ticklabels = {expst.stimdata.p2, round(expst.stimdata.p1, 3)};
        % histogram bins for each dimension
        expst.viz.binedges = cell(1, 2);
        expst.viz.binedges{1} = interp1(1:2:numel(expst.dimvals{1})*2, ...
                                      expst.dimvals{1}, ...
                                      0:2:numel(expst.dimvals{1})*2, ...
                                      "linear", "extrap");
        expst.viz.binedges{2} = interp1(1:2:numel(expst.dimvals{2})*2, ...
                                      log10(expst.dimvals{2}), ...
                                      0:2:numel(expst.dimvals{2})*2, ...
                                      "linear", "extrap");
    case "battery1"
        expst.stimlen = 16;
        %%%% assign stimulus names and values
        expst.stimabbrv = "B1";
        expst.dimnames = ["Spatial Frequency (cpd)", sprintf("Direction (%s)", char(176))];
        expst.dimabbrv = ["SF", "Dir"];
        expst.dimvals = {params.p2', params.p1'};
        expst.dimscales = ["log", "polar"];
        %%%% function handles for estimating peak parameters
        expst.calc.estfunc = {@sf_est, @dir_est};
        %%%% add visualization parameters
        expst.viz.window = [-30, 30];
        expst.viz.labels = {"Spatial Frequency (cpd)", sprintf("Direction (%s)", char(176))};
        expst.viz.ticks.cont = {log10(expst.stimdata.p2), expst.stimdata.p1};
        expst.viz.ticks.disc = {1:expst.stimdata.np2, 1:expst.stimdata.np1};
        expst.viz.ticklabels = {round(expst.stimdata.p2, 2), expst.stimdata.p1};
        % histogram bins for each dimension
        expst.viz.binedges = cell(1, 2);
        expst.viz.binedges{1} = interp1(1:2:numel(expst.dimvals{1})*2, ...
                                      log10(expst.dimvals{1}), ...
                                      0:2:numel(expst.dimvals{1})*2, ...
                                      "linear", "extrap");
        expst.viz.binedges{2} = interp1(1:2:numel(expst.dimvals{2})*2, ...
                                      expst.dimvals{2}, ...
                                      0:2:numel(expst.dimvals{2})*2, ...
                                      "linear", "extrap");
    case "battery2"
        expst.stimlen = 16;
        %%%% convert temporal period in frames to Hz
        params.p2 = 60 ./ params.p2;
        expst.stimdata = params;
        %%%% assign stimulus names and values
        expst.stimabbrv = "B2";
        expst.dimnames = ["Temporal Frequency (Hz)", sprintf("Direction (%s)", char(176))];
        expst.dimabbrv = ["TF", "Dir"];
        expst.dimvals = {flip(params.p2)', params.p1'};
        expst.dimscales = ["log", "polar"];
        %%%% function handles for estimating peak parameters
        expst.calc.estfunc = {@tf_est, @dir_est};
        %%%% add visualization parameters
        expst.viz.window = [-45, 45];
        expst.viz.labels = {"Temporal Frequency (Hz)", sprintf("Direction (%s)", char(176))};
        expst.viz.ticks.cont = {flip(log10(expst.stimdata.p2)), expst.stimdata.p1};
        expst.viz.ticks.cont{1}(2) = log10(60/76); %%%% move the 1 Hz tick
        expst.viz.ticks.disc = {1:expst.stimdata.np2, 1:expst.stimdata.np1};
        expst.viz.ticklabels = {flip(round(expst.stimdata.p2, 1)), expst.stimdata.p1};
        % histogram bins for each dimension
        expst.viz.binedges = cell(1, 2);
        tfvals = 60./[120, 76, 48, 31, 20, 12, 8];
        expst.viz.binedges{1} = 0.01 + interp1(1:2:numel(tfvals)*2, ...
                                       log10(tfvals), ...
                                       0:2:numel(tfvals)*2, ...
                                       "linear", "extrap");
        expst.viz.binedges{2} = interp1(1:2:numel(expst.dimvals{2})*2, ...
                                      expst.dimvals{2}, ...
                                      0:2:numel(expst.dimvals{2})*2, ...
                                      "linear", "extrap");
    case "battery3"
        expst.stimlen = 31;
        %%%% convert contrast in fraction to percent
        params.p2 = 100 .* params.p2;
        %%%% convert size from radius to diameter
        params.p3 = 2 .* params.p3;
        expst.stimdata = params;
        %%%% assign stimulus names and values
        expst.stimabbrv = "B3";
        expst.dimnames = ["Contrast (%)", sprintf("Direction (%s)", char(176)), sprintf("Size (%s)", char(176))];
        expst.dimabbrv = ["Ctrst", "Dir", "Size"];
        expst.dimvals = {params.p2', params.p1', params.p3'};
        expst.dimscales = ["log", "polar", "linear"];
        %%%% function handles for estimating peak parameters
        expst.calc.estfunc = {@contrast_est, @dir_est, @size_est};
        %%%% add visualization parameters
        expst.viz.window = [-30, 30];
        expst.viz.labels = {"Contrast (%)", sprintf("Direction (%s)", char(176)), sprintf("Size (%s)", char(176))};
        expst.viz.ticks.cont = {log10(expst.stimdata.p2), expst.stimdata.p1, expst.stimdata.p3};
        expst.viz.ticks.cont{1}(1) = expst.viz.ticks.cont{1}(2) - mean(diff(expst.viz.ticks.cont{1}(2:end)));
        expst.viz.ticks.cont{3}(1) = 0; %%%% change the location of the first size tick to center it
        expst.viz.ticks.disc = {1:expst.stimdata.np2, 1:expst.stimdata.np1, 1:expst.stimdata.np3};
        expst.viz.ticklabels = {round(expst.stimdata.p2, 0), expst.stimdata.p1, expst.stimdata.p3};
        expst.viz.binedges{1} = interp1(1:2:numel(expst.dimvals{1})*2, ...
                                      log10(expst.dimvals{1}), ...
                                      0:2:numel(expst.dimvals{1})*2, ...
                                      "linear", "extrap");
        % account for log10 of zero
        expst.viz.binedges{1}(2) = expst.viz.binedges{1}(3) - mean(diff(expst.viz.binedges{1}(3:end)));
        expst.viz.binedges{1}(1) = expst.viz.binedges{1}(3) - 2 * mean(diff(expst.viz.binedges{1}(3:end)));
        expst.viz.binedges{2} = interp1(1:2:numel(expst.dimvals{2})*2, ...
                                      expst.dimvals{2}, ...
                                      0:2:numel(expst.dimvals{2})*2, ...
                                      "linear", "extrap");
        expst.viz.binedges{3} = -4.99:10:65.01; %%%% make size bins equal
        % interp1(1:2:numel(expst.dimvals{3})*2, ...
        %                               expst.dimvals{3}, ...
        %                               0:2:numel(expst.dimvals{3})*2, ...
        %                               "linear", "extrap");
    case "battery4"
        expst.stimlen = 16;
        %%%% convert contrast in fraction to percent and multiply by 2
        params.p3 = 100 .* params.p3 .* 2;
        params.p4 = 100 .* params.p4 .* 2;
        %%%% convert size from radius to diameter
        params.p5 = 2 .* params.p5;
        expst.stimdata = params;
        %%%% find center, cross, iso, blank stimuli
        stimparams = params.param(params.idx, :);
        center_blank = stimparams(:, 3) == 0;
        surround_blank = stimparams(:, 4) == 0;
        center_up = stimparams(:, 1) == 0; center_down = stimparams(:, 1) == 180;
        center_left = stimparams(:, 1) == 90; center_right = stimparams(:, 1) == 270;
        surround_up = stimparams(:, 2) == 0; surround_down = stimparams(:, 2) == 180;
        surround_left = stimparams(:, 2) == 90; surround_right = stimparams(:, 2) == 270;
        center_ud = ismember(stimparams(:, 1), [0, 180]);
        center_lr = ismember(stimparams(:, 1), [90, 270]);
        surround_ud = ismember(stimparams(:, 2), [0, 180]);
        surround_lr = ismember(stimparams(:, 2), [90, 270]);
        expst.stimdata.context.center = ~center_blank & surround_blank;
        expst.stimdata.context.surround = center_blank & ~surround_blank;
        expst.stimdata.context.cross = ((center_ud & surround_lr) | (center_lr & surround_ud)) & ...
                                        (~center_blank & ~surround_blank);
        expst.stimdata.context.iso = ((center_up & surround_up) | (center_down & surround_down) | ...
                                      (center_left & surround_left) | (center_right & surround_right)) & ...
                                      (~center_blank & ~surround_blank);
        expst.stimdata.context.opp = ((center_up & surround_down) | (center_down & surround_up) | ...
                                      (center_left & surround_right) | (center_right & surround_left)) & ...
                                      (~center_blank & ~surround_blank);
        expst.stimdata.context.blank = center_blank & surround_blank;
        %%%% assign stimulus names and values
        expst.stimabbrv = "B4";
        expst.dimnames = [sprintf("Center Direction (%s)", char(176)), sprintf("Surround Direction (%s)", char(176)), ...
                        "Center Contrast (%)", "Surround Contrast (%)", sprintf("Size (%s)", char(176))];
        expst.dimabbrv = ["Dir_Cent", "Dir_Surr", "Ctrst_Cent", "Ctrst_Surr", "Size"];
        expst.dimvals = {params.p2', params.p1', params.p3', params.p4', params.p5'};
        expst.dimscales = ["polar", "polar", "log", "log", "linear"];
        %%%% function handles for estimating peak parameters
        expst.calc.estfunc = {@dir_est, @dir_est, @contrast_est, @contrast_est, @size_est};
        %%%% add visualization parameters
        expst.viz.window = [-30, 30];
        expst.viz.labels = {sprintf("Center Direction (%s)", char(176)), sprintf("Surround Direction (%s)", char(176)), ...
                          "Center Contrast (%)", "Surround Contrast (%)", sprintf("Size (%s)", char(176))};
        expst.viz.ticks.cont = {expst.stimdata.p2, expst.stimdata.p1, log10(expst.stimdata.p3), log10(expst.stimdata.p4), expst.stimdata.p5};
        expst.viz.ticks.disc = {1:expst.stimdata.np2, 1:expst.stimdata.np1, 1:expst.stimdata.np3, 1:expst.stimdata.np4, 1:expst.stimdata.np5};
        expst.viz.ticklabels = {expst.stimdata.p2, expst.stimdata.p1, expst.stimdata.p3, expst.stimdata.p4, expst.stimdata.p5};
        expst.viz.binedges{1} = interp1(1:2:numel(expst.dimvals{1})*2, ...
                                        expst.dimvals{1}, ...
                                        0:2:numel(expst.dimvals{1})*2, ...
                                        "linear", "extrap");
        expst.viz.binedges{2} = interp1(1:2:numel(expst.dimvals{2})*2, ...
                                        expst.dimvals{2}, ...
                                        0:2:numel(expst.dimvals{2})*2, ...
                                        "linear", "extrap");
        expst.viz.binedges{3} = interp1(1:2:numel(expst.dimvals{3})*2, ...
                                        log10(expst.dimvals{3}), ...
                                        0:2:numel(expst.dimvals{3})*2, ...
                                        "linear", "extrap");
        % account for log10 of zero
        expst.viz.binedges{3}(2) = expst.viz.binedges{3}(3) - mean(diff(expst.viz.binedges{3}(3:end)));
        expst.viz.binedges{3}(1) = expst.viz.binedges{3}(3) - 2 * mean(diff(expst.viz.binedges{3}(3:end)));
        expst.viz.binedges{4} = interp1(1:2:numel(expst.dimvals{4})*2, ...
                                        log10(expst.dimvals{4}), ...
                                        0:2:numel(expst.dimvals{4})*2, ...
                                        "linear", "extrap");
        % account for log10 of zero
        expst.viz.binedges{4}(2) = expst.viz.binedges{4}(3) - mean(diff(expst.viz.binedges{4}(3:end)));
        expst.viz.binedges{4}(1) = expst.viz.binedges{4}(3) - 2 * mean(diff(expst.viz.binedges{4}(3:end)));
        expst.viz.binedges{5} = interp1(1:2:numel(expst.dimvals{5})*2, ...
                                        expst.dimvals{5}, ...
                                        0:2:numel(expst.dimvals{5})*2, ...
                                        "linear", "extrap");
end

%%%% get different forms of calcium data
traces = num2cell(traces, 1);
spikes = num2cell(spikes, 1);
% add dF/F based on stimuli
traces_mat = horzcat(traces{:});
spikes_mat = horzcat(spikes{:});
baseline = zeros(expst.nstims, expst.ncells);
for s = 1:expst.nstims
    baseline(s, :) = mean(traces_mat(frameon(s) : frameon(s) + 5, :), 1);
end
if stim_type == "randorisf", baseline = traces_mat; end
dff = (traces_mat - mean(baseline, 1)) ./ mean(baseline, 1);
dFF = num2cell(dff, 1);

%%%% add stims
if stim_type == "battery4"
    param = params.param(params.idx, :) .* [1, 1, 200, 200, 1]; % map the contrasts to percentages
else
    param = params.param(params.idx, :);
end
onsets = num2cell(frameon);
frameoff = frameon + expst.stimlen;
offsets = num2cell(frameoff);
lens = num2cell(frameoff - frameon);
idxs = num2cell(params.idx');

expst.stims = struct("param", num2cell(param, 2), ...
                     "onset", onsets, ...
                     "offset", offsets, ...
                     "nframes", lens, ...
                     "idx", idxs);

%%%% assign stimuli as having significant (5 pixels) eye divergence
eye_pos = vertcat(eye.Centroid);
eye_pos = eye_pos(1:expst.nframes, :);
eye_pos_logical = (eye_pos(:, 1) > eye_center(1) - eye_radius) & ...
                  (eye_pos(:, 1) < eye_center(1) + eye_radius) & ...
                  (eye_pos(:, 2) > eye_center(2) - eye_radius) & ...
                  (eye_pos(:, 2) < eye_center(2) + eye_radius);
% assign stimuli as running or not running
stimeye = cell(expst.nstims, 1);
for st = 1:expst.nstims
    onset = expst.stims(st).onset;
    offset = expst.stims(st).offset;
    if all(eye_pos_logical(onset:offset))
        stimeye{st} = true;
    else
        stimeye{st} = false;
    end
end
[expst.stims.eye] = stimeye{:};

%%%% iterate over different forms of data and compute
for data_form = 1:2

    if data_form == 1,     data = spikes_mat; expst.dtype = "spikes"; window = [5, 0];
    elseif data_form == 2, data = dff;        expst.dtype = "dFF";    window = [8, 8];
    end
    
    full_grid = repelem(reshape(1:9, [3, 3])', ops.Ly / 3, ops.Lx / 3);
    xy = zeros(expst.ncells, 2);
    grid = zeros(expst.ncells, 1);
    for c = 1:expst.ncells
        loc = stat{c}.med;
        xy(c, :) = loc;
        grid(c) = full_grid(loc(1), loc(2));
    end
 
    expst.cells = struct("data", num2cell(data, 1), ...
                         "xy", num2cell(xy, 2)', ...
                         "grid", num2cell(grid, 2)');

    %%%% add kernels
    if stim_type ~= "randorisf"
        kernel_out = compute_kernel(data, frameon, frameoff, params, window, vertcat(expst.stims.eye), 0:4);
    else
        [~, kernel_out, stats] = compute_kernel_orisf(data, frameon, frameoff, params);
    end
    if stim_type == "battery2", kernel_out = flip(kernel_out, 1); end
    kerns = permute(num2cell(kernel_out, 1:expst.ndims), flip(1:expst.ndims+1));
    peaks = squeeze(num2cell(max(kernel_out, [], 1:expst.ndims, "omitnan")));
    valleys = squeeze(num2cell(min(kernel_out, [], 1:expst.ndims, "omitnan")));
    % find tuning curves at max for each dimension
    curves = cell(numel(kerns), expst.ndims);
    peakidxs = cell(numel(kerns), 1);
    peakvals = cell(numel(kerns), 1);
    estvals = cell(numel(kerns), expst.ndims);
    for k = 1:numel(kerns)
        idx = nan(1, expst.ndims);
        vals = nan(1, expst.ndims);
        for d = 1:expst.ndims
            kernshift = permute(kerns{k}, [d, setdiff(1:expst.ndims, d)]);
            [idx(d), ~] = find(kernshift == peaks{k}, 1);
        end
        peakidxs{k} = idx;
        for v = 1:expst.ndims
            vals(v) = expst.dimvals{v}(idx(v));
        end
        peakvals{k} = vals;
        for d = 1:expst.ndims
            kernshift = permute(kerns{k}, [d, setdiff(1:expst.ndims, d)]);
            restdims = idx(setdiff(1:expst.ndims, d));
            switch expst.ndims
                case 2
                    curves{k, d} = reshape(kernshift(:, restdims(1)), 1, []);
                case 3
                    curves{k, d} = reshape(kernshift(:, restdims(1), restdims(2)), 1, []);
                case 4
                    curves{k, d} = reshape(kernshift(:, restdims(1), restdims(2), restdims(3)), 1, []);
                case 5
                    curves{k, d} = reshape(kernshift(:, restdims(1), restdims(2), restdims(3), restdims(4)), 1, []);
            end
            estvals{k, d} = expst.calc.estfunc{d}(curves{k, d}, expst.dimvals{d});
        end
    end

    expst.kernels = struct("kern", kerns', ...
                           "kernsmooth", cellfun(@(x) smooth_kernel(x, 3), kerns, 'UniformOutput', false)', ...
                           "curve", num2cell(curves, 2)', ...
                           "peakidx", peakidxs', ...
                           "peakval", peakvals', ...
                           "estval", num2cell(cell2mat(estvals), 2)', ...
                           "peak", peaks', ...
                           "valley", valleys');
    
    %%%% add PSTHs
    if stim_type ~= "randorisf"
        [psth_kern, psth_all] = compute_psth(data, frameon, frameoff, params, expst.viz.window, 0:4);
        psthmats = squeeze(num2cell(psth_all, [1, 2]));
        psthkerns = squeeze(num2cell(psth_kern, 1:expst.ndims));
        peaks = cellfun(@(x) max(x, [], "all"), psthmats, "UniformOutput", false);
        valleys = cellfun(@(x) min(x, [], "all"), psthmats, "UniformOutput", false);
        windows = repmat({expst.viz.window}, [expst.ncells, 1]);
    
        expst.psths = struct("mat", psthmats', ...
                             "kern", psthkerns', ...
                             "peak", peaks', ...
                             "valley", valleys', ...
                             "window", windows');
    end

%%%% add stimulus-non-specific statistics
% skewness
expst.stats = struct("skewness", num2cell(skewness(data), 1)');
% pupil diameter
eye_data = 2 .* vertcat(eye.Radius); % use diameter
expst.eye.data = medfilt1(eye_data(1:expst.nframes), 15, "truncate"); % median filter to remove blinking artifacts
expst.eye.pos = eye_pos;
expst.eye.logical = logical(eye_pos_logical);
eye_corr = num2cell(corr(horzcat(expst.cells.data), expst.eye.data));
[expst.stats.eye_corr] = eye_corr{:};
% running (absolute first derivative of quadrature)
quad_data = double(quad');
expst.quad.data = movmean(quad_data, 15);
expst.quad.logical = expst.quad.data > 0.5 * std(expst.quad.data);
% assign stimuli as running or not running
stimrun = cell(expst.nstims, 1);
for st = 1:expst.nstims
    onset = expst.stims(st).onset;
    offset = expst.stims(st).offset;
    if all(expst.quad.logical(onset:offset))
        stimrun{st} = 2;
    elseif any(expst.quad.logical(onset:offset))
        stimrun{st} = 1;
    else
        stimrun{st} = 0;
    end
end
[expst.stims.run] = stimrun{:};
% compute full-session correlation with running
quad_corr = num2cell(corr(horzcat(expst.cells.data), expst.quad.data));
[expst.stats.quad_corr] = quad_corr{:};

%%%% calculate stimulus-specific statistics
switch stim_type
    case "randorisf"
        % assign stats computed in compute_kernel_orisf
        [expst.stats.SNR] = stats.SNR;
        % pause;
        [expst.stats.tmax] = stats.tmax;
        [expst.stats.F1F0] = stats.F1F0;
        % compute the population SNR threshold
        expst.calc.SNR_thr = SNR_thr(vertcat(expst.stats.SNR), vertcat(expst.stats.tmax));
        % compute circular variance
        curves = vertcat(expst.kernels.curve);
        CV = cellfun(@(x) circ_var(x, expst.dimvals{1}), curves(:, 1), "UniformOutput", false);
        [expst.stats.CV] = CV{:};
    case "battery1"
        curves = vertcat(expst.kernels.curve);
        DSI = cellfun(@(x) dir_slc_idx(x, expst.dimvals{2}), curves(:, 2), "UniformOutput", false);
        [expst.stats.DSI] = DSI{:};
    case "battery2"
        curves = vertcat(expst.kernels.curve);
        DSI = cellfun(@(x) dir_slc_idx(x, expst.dimvals{2}), curves(:, 2), "UniformOutput", false);
        [expst.stats.DSI] = DSI{:};
    case "battery3"
        curves = vertcat(expst.kernels.curve);
        DSI = cellfun(@(x) dir_slc_idx(x, expst.dimvals{2}), curves(:, 2), "UniformOutput", false);
        [expst.stats.DSI] = DSI{:};
        % stationary/running size tuning
        stationary = cell(expst.ncells, 1); running = cell(expst.ncells, 1);
        PSTHs = {expst.psths.mat};
        for c = 1:numel(PSTHs)
        % get stimuli where eye is looking in the right spot
            sub_idx = vertcat(expst.stims.eye);
            [s, r] = size_est_quad(PSTHs{c}, expst.stimdata, vertcat(expst.stims.run), expst.viz.window, [0, 4], [10, 25], sub_idx);
            stationary{c} = s; running{c} = r;
        end
        [expst.stats.size_stationary] = stationary{:};
        [expst.stats.size_running] = running{:};
    case "battery4"
        % define kernel statistics
        kern_max = @(k) num2cell(reshape(max(k, [], 1:4, "omitnan"), [], expst.ncells)', 2);
        kern_max_all = @(k) num2cell(repmat(reshape(max(k, [], 1:5, "omitnan"), [], expst.ncells)', [1, size(k, 5)]), 2);
        kern_mean = @(k) num2cell(reshape(mean(k, 1:4, "omitnan"), [], expst.ncells)', 2);
        kern_mean_all = @(k) num2cell(repmat(reshape(mean(k, 1:5, "omitnan"), [], expst.ncells)', [1, size(k, 5)]), 2);
        kern_std = @(k) num2cell(reshape(std(k, 0, 1:4, "omitnan"), [], expst.ncells)', 2);
        kern_std_all = @(k) num2cell(repmat(reshape(std(k, 0, 1:5, "omitnan"), [], expst.ncells)', [1, size(k, 5)]), 2);
        mod_idx = @(c1, c2) num2cell((cell2mat(c1) - cell2mat(c2)) ./ (cell2mat(c1) + cell2mat(c2)), 2);
        iti = @(c, s, i) num2cell(((cell2mat(c) - cell2mat(s)) ./ 2 .* ((cell2mat(s) - cell2mat(i)) + (cell2mat(c) - cell2mat(i)))) + 0.5, 2);
        % center response
        kernel_center = compute_kernel(data, frameon, frameoff, params, window, expst.stimdata.context.center & vertcat(expst.stims.eye), 0:4);
        center_max = kern_max(kernel_center); center_mean = kern_mean(kernel_center); center_std = kern_std(kernel_center);
        % surround response
        kernel_surround = compute_kernel(data, frameon, frameoff, params, window, expst.stimdata.context.surround & vertcat(expst.stims.eye), 0:4);
        surround_max = kern_max(kernel_surround); surround_mean = kern_mean(kernel_surround); surround_std = kern_std(kernel_surround);
        % iso response
        kernel_iso = compute_kernel(data, frameon, frameoff, params, window, expst.stimdata.context.iso & vertcat(expst.stims.eye), 0:4);
        iso_max = kern_max(kernel_iso); iso_mean = kern_mean(kernel_iso); iso_std = kern_std(kernel_iso);
        iso_all_max = kern_max_all(kernel_iso); iso_all_mean = kern_mean_all(kernel_iso); iso_all_std = kern_std_all(kernel_iso); % average across size too
        % cross response
        kernel_cross = compute_kernel(data, frameon, frameoff, params, window, expst.stimdata.context.cross & vertcat(expst.stims.eye), 0:4);
        cross_max = kern_max(kernel_cross); cross_mean = kern_mean(kernel_cross); cross_std = kern_std(kernel_cross);
        % blank response
        kernel_blank = compute_kernel(data, frameon, frameoff, params, window, expst.stimdata.context.blank & vertcat(expst.stims.eye), 0:4);
        blank_max = kern_max(kernel_blank); blank_mean = kern_mean(kernel_blank); blank_std = kern_std(kernel_blank);
        % modulation index (cross / iso)
        cri_mod_max = mod_idx(cross_max, iso_all_max);
        cri_mod_mean = mod_idx(cross_mean, iso_all_mean);
        % modulation index (center / surround)
        cs_mod_max = mod_idx(cross_max, surround_max);
        cs_mod_mean = mod_idx(cross_mean, surround_mean);
        % modulation index (center / iso)
        ci_mod_max = mod_idx(center_max, iso_max);
        ci_mod_mean = mod_idx(center_mean, iso_mean);
        % ITI
        ITI = iti(center_max, surround_max, iso_all_max);
        % center/surround size tuning curve correlation
        cs_size_corr = cellfun(@(x, y) corrcoef(x, y), center_mean, surround_mean, "UniformOutput", false);
        cs_size_corr = cellfun(@(x) x(2, 1), cs_size_corr, "UniformOutput", false);
        %%%% assign to stats field
        [expst.stats.center_max] = center_max{:};
        [expst.stats.center_mean] = center_mean{:};
        [expst.stats.center_std] = center_std{:};
        [expst.stats.surround_max] = surround_max{:};
        [expst.stats.surround_mean] = surround_mean{:};
        [expst.stats.surround_std] = surround_std{:};
        [expst.stats.iso_max] = iso_max{:};
        [expst.stats.iso_all_max] = iso_all_max{:};
        [expst.stats.iso_mean] = iso_mean{:};
        [expst.stats.iso_all_mean] = iso_all_mean{:};
        [expst.stats.iso_std] = iso_std{:};
        [expst.stats.iso_all_std] = iso_all_std{:};
        [expst.stats.cross_max] = cross_max{:};
        [expst.stats.cross_mean] = cross_mean{:};
        [expst.stats.cross_std] = cross_std{:};
        [expst.stats.blank_max] = blank_max{:};
        [expst.stats.blank_mean] = blank_mean{:};
        [expst.stats.blank_std] = blank_std{:};
        [expst.stats.cri_mod_max] = cri_mod_max{:};
        [expst.stats.cri_mod_mean] = cri_mod_mean{:};
        [expst.stats.cs_mod_max] = cs_mod_max{:};
        [expst.stats.cs_mod_mean] = cs_mod_mean{:};
        [expst.stats.ci_mod_max] = ci_mod_max{:};
        [expst.stats.ci_mod_mean] = ci_mod_mean{:};
        [expst.stats.ITI] = ITI{:};
        [expst.stats.cs_size_corr] = cs_size_corr{:};
end

%%%% generate tables for all stimuli (estimated tuning, running, pupil, calcium dynamics)
% add tuning parameter estimates
estvals = vertcat(expst.kernels.estval);
expst.table = array2table(estvals, VariableNames=strcat(expst.dimabbrv + "_" + expst.stimabbrv + "_" + expst.dtype));
% add table property for storing feature scale (linear/log)
expst.table = addprop(expst.table, "FeatureScale", "variable");
expst.table.Properties.CustomProperties.FeatureScale = expst.dimscales;
% add table property for storing bin edges for feature histograms
expst.table = addprop(expst.table, "BinEdges", "variable");
expst.table.Properties.CustomProperties.BinEdges = expst.viz.binedges;
% add table property for storing ticks for feature plots
expst.table = addprop(expst.table, "Ticks", "variable");
expst.table.Properties.CustomProperties.Ticks = expst.viz.ticks.cont;
% add table property for storing ticklabels for feature plots
expst.table = addprop(expst.table, "TickLabels", "variable");
expst.table.Properties.CustomProperties.TickLabels = expst.viz.ticklabels;
% add response magnitudes?
stats = ["skewness", "eye_corr", "quad_corr"];
binedges = {-2:10, -0.5:0.1:0.5, -0.5:0.1:0.5};
ticks = {0:2:10, -0.5:0.25:0.5, -0.5:0.25:0.5};
for s = 1:numel(stats)
    expst.table.(stats(s) + "_" + expst.stimabbrv + "_" + expst.dtype) = vertcat(expst.stats.(stats(s)));
    expst.table.Properties.CustomProperties.FeatureScale{end} = 'linear';
    expst.table.Properties.CustomProperties.BinEdges{end} = binedges{s};
    expst.table.Properties.CustomProperties.Ticks{end} = ticks{s};
    expst.table.Properties.CustomProperties.TickLabels{end} = ticks{s};
end
%%%% generate tables for each stimulus
switch stim_type
    case "randorisf"
        stats = ["SNR", "tmax", "F1F0", "CV"];
        binedges = {0:0.5:6, 1:20, 0:0.2:2, 0:0.1:1};
        ticks = {0:6, 1:4:20, 0:0.5:2, 0:0.25:1};
        for s = 1:numel(stats)
            expst.table.(stats(s) + "_" + expst.stimabbrv + "_" + expst.dtype) = vertcat(expst.stats.(stats(s)));
            expst.table.Properties.CustomProperties.FeatureScale{end} = 'linear';
            expst.table.Properties.CustomProperties.BinEdges{end} = binedges{s};
            expst.table.Properties.CustomProperties.Ticks{end} = ticks{s};
            expst.table.Properties.CustomProperties.TickLabels{end} = ticks{s};
        end
    case "battery1"
        % direction tuning, spatial frequency tuning added above
        stats = ["DSI"];
        binedges = {0:0.1:1.1};
        ticks = {0:0.2:1};
        for s = 1:numel(stats)
            expst.table.(stats(s) + "_" + expst.stimabbrv + "_" + expst.dtype) = vertcat(expst.stats.(stats(s)));
            expst.table.Properties.CustomProperties.FeatureScale{end} = 'linear';
            expst.table.Properties.CustomProperties.BinEdges{end} = binedges{s};
            expst.table.Properties.CustomProperties.Ticks{end} = ticks{s};
            expst.table.Properties.CustomProperties.TickLabels{end} = ticks{s};
        end
    case "battery2"
        % direction tuning, temporal frequency tuning added above
        stats = ["DSI"];
        binedges = {0:0.1:1.1};
        ticks = {0:0.2:1};
        for s = 1:numel(stats)
            expst.table.(stats(s) + "_" + expst.stimabbrv + "_" + expst.dtype) = vertcat(expst.stats.(stats(s)));
            expst.table.Properties.CustomProperties.FeatureScale{end} = 'linear';
            expst.table.Properties.CustomProperties.BinEdges{end} = binedges{s};
            expst.table.Properties.CustomProperties.Ticks{end} = ticks{s};
            expst.table.Properties.CustomProperties.TickLabels{end} = ticks{s};
        end
    case "battery3"
        % direction tuning, contrast tuning, size tuning added above
        stats = ["DSI"];
        binedges = {0:0.1:1.1};
        ticks = {0:0.2:1};
        for s = 1:numel(stats)
            expst.table.(stats(s) + "_" + expst.stimabbrv + "_" + expst.dtype) = vertcat(expst.stats.(stats(s)));
            expst.table.Properties.CustomProperties.FeatureScale{end} = 'linear';
            expst.table.Properties.CustomProperties.BinEdges{end} = binedges{s};
            expst.table.Properties.CustomProperties.Ticks{end} = ticks{s};
            expst.table.Properties.CustomProperties.TickLabels{end} = ticks{s};
        end
    case "battery4"
        % center response, surround response, iso response, cross response,
        % modulation index (at preferred and all stimuli)
        stats = ["cri_mod_max", "cri_mod_mean", "cs_mod_max", "cs_mod_mean", "ci_mod_max", "ci_mod_mean", "ITI", "cs_size_corr"];
        binedges = {0:0.1:1, -1:0.2:1, 0:0.1:1, -1:0.2:1, 0:0.1:1, -1:0.2:1, 0:0.1:1, -1:0.2:1};
        ticks = {0:0.25:1, -1:0.5:1, 0:0.25:1, -1:0.5:1, 0:0.25:1, -1:0.5:1, 0:0.5:1, -1:0.5:1};
        for s = 1:numel(stats)
            st = vertcat(expst.stats.(stats(s)));
            if size(st, 2) > 1
                expst.table.(stats(s) + "_" + expst.stimabbrv + "_" + expst.dtype) = st(:, 3); %%%% grab 20 degrees for table
            else
                expst.table.(stats(s) + "_" + expst.stimabbrv + "_" + expst.dtype) = st;
            end
            expst.table.Properties.CustomProperties.FeatureScale{end} = 'linear';
            expst.table.Properties.CustomProperties.BinEdges{end} = binedges{s};
            expst.table.Properties.CustomProperties.Ticks{end} = ticks{s};
            expst.table.Properties.CustomProperties.TickLabels{end} = ticks{s};
        end
end

%%%% assign structure to the appropriate form of the data
if data_form == 1,     spikes = expst;
elseif data_form == 2, dFF = expst;
end

end

    function kernel_smoothed = smooth_kernel(kernel, nbins)
        if ismatrix(kernel)
            h = fspecial('gauss', nbins, 1);
            % extend kernel by replicating edges for filtering
            kernel2 = [kernel(1, :); kernel; kernel(end, :)];
            kernel3 = [kernel2(:, 1), kernel2, kernel2(:, end)];
            kernel_smoothed = filter2(h, kernel3, 'valid');
        else
            h = fspecial3('gauss', [nbins, nbins, nbins]);
            kernel_smoothed = imfilter(kernel, h, 'replicate');
        end
    end

    function est = ori_est(curve, dimvals)
        est = rad2deg(angle(sum(curve .* exp(1i * dimvals * 2 * pi/180))) / 2);
        if est < 0
            est = est + 180;
        end
    end

    function cv = circ_var(curve, dimvals)
        cv = 1 - abs(sum(curve .* exp(1i * 2 * pi * dimvals/180)) ./ sum(curve));
    end

    function est = sf_est(curve, dimvals)
        resp_sf = curve - max(curve) .* .75; % clip tails...
        resp_sf(resp_sf < 0) = 0;
        est = 10 ^ (sum(resp_sf .* log10(dimvals)) / sum(resp_sf));
    end

    function thr = SNR_thr(SNR, tmax)
        tmax_idx = tmax >= 5 & tmax <= 11;
        SNR_mean = mean(SNR(~tmax_idx), "omitnan");
        SNR_std = std(SNR(~tmax_idx), "omitnan");
        thr = SNR_mean + 3 * SNR_std;
        if isnan(thr), thr = 1.5; end
    end

    function est = dir_est(curve, dimvals)
        curve(isnan(curve)) = 0;
        est = rad2deg(angle(sum(curve .* exp(1i * dimvals * pi/180))));
        if est < 0
            est = est + 360;
        end
    end

    function DSI = dir_slc_idx(curve, dimvals)
        peak = find(curve == max(curve), 1);
        Rp = curve(peak);
        if dimvals(peak) < 180
            Rn = curve(dimvals == (dimvals(peak) + 180));
        else
            Rn = curve(dimvals == (dimvals(peak) - 180));
        end
        DSI = (Rp - Rn) / (Rp + Rn);
    end

    function est = tf_est(curve, dimvals)
        resp_tf = curve - max(curve) .* .75; % clip tails...
        resp_tf(resp_tf < 0) = 0;
        est = 10 ^ (sum(resp_tf .* log10(dimvals)) / sum(resp_tf));
    end

    function est = contrast_est(curve, dimvals)
        resp_contrast = curve - max(curve) .* .75; % clip tails...
        resp_contrast(resp_contrast < 0) = 0;
        est = 10 ^ (sum(resp_contrast .* log10(dimvals + 1)) / sum(resp_contrast));
        est = est - 1;
    end

    function est = size_est(curve, dimvals)
        resp_size = curve - max(curve) .* .75; % clip tails...
        resp_size(resp_size < 0) = 0;
        est = sum(resp_size .* dimvals) / sum(resp_size);
    end

    function [stationary, running] = size_est_quad(psth, stimdata, run, window, baseline, response, sub_idx)
        if isnan(sub_idx), sub_idx = true(numel(run), 1); end
        wdw = window(1):window(2);
        trial_baseline = mean(psth(sub_idx, ismember(wdw, baseline(1):baseline(2))), 2);
        trial_response = mean(psth(sub_idx, ismember(wdw, response(1):response(2))), 2);
        trial_subtract = trial_response - trial_baseline;
        size_avg = zeros(stimdata.np3, 1);
        sizes = stimdata.param(stimdata.idx(sub_idx), 3);
        for sz = 1:stimdata.np3
            size_avg(sz) = mean(trial_subtract(sizes == stimdata.p3(sz)));
        end
        stationary = zeros(1, stimdata.np3); running = zeros(1, stimdata.np3);
        stationary_count = zeros(1, stimdata.np3); running_count = zeros(1, stimdata.np3);
        run = run(sub_idx);
        stimidx = stimdata.stimidx(sub_idx, :);
        for stm = 1:size(stimidx, 1)
            trial_data = trial_subtract(stm);
            if run(stm) == 0 % not running
                stationary(stimidx(stm, 3)) = stationary(stimidx(stm, 3)) + trial_data;
                stationary_count(stimidx(stm, 3)) = stationary_count(stimidx(stm, 3)) + 1;
            elseif run(stm) > 1 % running the whole time
                running(stimidx(stm, 3)) = running(stimidx(stm, 3)) + trial_data;
                running_count(stimidx(stm, 3)) = running_count(stimidx(stm, 3)) + 1;
            end
        end
        stationary = stationary ./ stationary_count;
        running = running ./ running_count;
    end

end