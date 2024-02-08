function navg = compute_kernel(pop, frameon, frameoff, params, window, sub_idx, subtract_fr)

if nargin < 5, window = [0, 0]; end
if nargin < 6, sub_idx = NaN; end
if nargin < 7, subtract_fr = NaN; end

if ~isnan(sub_idx)
    frameon = frameon(sub_idx);
    frameoff = frameoff(sub_idx);
    params.idx = params.idx(sub_idx);
    params.stimidx = params.stimidx(sub_idx, :);
end

nstim = numel(frameon);
ncells = size(pop, 2);

vmean = zeros(nstim, ncells);
for i = 1:nstim                                     % for each stim index
    wdw = frameon(i) + window(1) : frameoff(i) + window(2);
    v = mean(pop(wdw, :));
    if ~isnan(subtract_fr)
        baseline = frameon(i) + subtract_fr;
        v = v - mean(pop(baseline, :));
    end
    vmean(i,:) = v;
end

idxs = unique(params.idx);
M = numel(idxs);
navg = nan([size(params.pp1) ncells]);             % and the average (across repeats) response to each
for k = 1:M
    j = params.idx == idxs(k);                     % find all the times we presented stim k
    s = unique(params.stimidx(j, :), "rows");      % make sure there's only one stim
    switch numel(s)
        case 2
            navg(s(2), s(1), :) = mean(vmean(j, :));                  % average the response
        case 3
            navg(s(2), s(1), s(3), :) = mean(vmean(j, :));                  % average the response
        case 4
            navg(s(2), s(1), s(3), s(4), :) = mean(vmean(j, :));                  % average the response
        case 5
            navg(s(2), s(1), s(3), s(4), s(5), :) = mean(vmean(j, :));                  % average the response
    end
end

end