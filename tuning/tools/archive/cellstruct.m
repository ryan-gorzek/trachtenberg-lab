function exp = expstruct(traces_mat, spikes_mat, frameon, frameoff, params, quad, eye)

exp = [];
exp.ncells = size(traces_mat, 2);
exp.nstims = numel(frameon);

% create cells
traces = cell(1, exp.ncells);
for t = 1:exp.ncells, traces{t} = traces_mat(:, t); end

spikes = cell(1, exp.ncells);
for s = 1:exp.ncells, spikes{s} = spikes_mat(:, s); end

exp.cells = struct("traces", traces, ...
                   "spikes", spikes);

% create stims
onsets = cell(1, exp.nstims);
for o = 1:exp.nstims, onsets{o} = frameon(o); end

offsets = cell(1, exp.nstims);
for o = 1:exp.nstims, offsets{o} = frameoff(o); end

exp.stims = struct("onsets", onsets, ...
                   "offsets", offsets);

end