function navg = compute_kernel(pop, frameon, frameoff, params, window)

if nargin < 5, window = [0, 0]; end

nstim = numel(frameon);
ncells = size(pop, 2);

vmean = zeros(nstim, ncells);
for i = 1:nstim                                     % for each stim index
    wdw = frameon(i) + window(1) : frameoff(i) + window(2);
    v = mean(pop(wdw, :));
    vmean(i,:) = v;
end

M = max(params.idx);
navg = zeros(M, ncells);                            % and the average (across repeats) response to each
for k = 1:M
    j = params.idx == k;                            % find all the times we presented stim k
    navg(k,:) = mean(vmean(j, :));                  % average the response
end

navg = reshape(navg, [size(params.pp1) ncells]);    % reshape so it has the right matrix form

end