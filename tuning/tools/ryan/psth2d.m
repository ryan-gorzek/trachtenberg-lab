function npsth = compute_psth(pop, frameon, frameoff, params, window)

if nargin < 5, window = [-15, 15]; end

nstim = numel(frameon);
ncells = size(pop, 2);
nframes = numel(max(frameon + window(2) - frameon + window(1)));

vmean = nan(nstim, nframes, ncells);
for i = 1:nstim                                        % for each stim index
    wdw = frameon(i) + window(1) : frameon(i) + window(2);
    v = pop(wdw, :);
    vmean(i, 1:numel(wdw), :) = v;
end

M = max(params.idx);
npsth = cell(M, ncells);                               % and the average (across repeats) response to each
for k = 1:M
    j = params.idx == k;                               % find all the times we presented stim k
    npsth(k, :) = num2cell(vmean(j, :, :), [1, 2]);    % average the response
end

npsth = reshape(npsth, [size(params.pp1) ncells]);     % reshape so it has the right matrix form

end