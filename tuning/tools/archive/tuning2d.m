function navg = tuning2d(pop, frameon, frameoff, param, idx, window)

if nargin < 6, window = [0, 0]; end

nstim = numel(frameon);
ncells = size(pop, 2);

vmean = zeros(nstim, ncells);
for i = 1:nstim                                         % for each stim index
    v = mean(pop(frameon(i) + window(1) : frameoff(i) + window(2), :));
    vmean(i,:) = v;
end

M = max(idx);
p1 = unique(param(:, 1)); %#ok<*USENS> 
p2 = unique(param(:, 2));
[pp1, ~] = meshgrid(p1,p2);
navg = zeros(M, ncells);                           % and the average (across repeats) response to each
for k = 1:M
    j = idx == k;                                       % find all the times we presented stim k
    navg(k,:) = mean(vmean(j, :));                  % average the response
end

navg = reshape(navg, [size(pp1) ncells]);     % reshape so it has the right matrix form

end