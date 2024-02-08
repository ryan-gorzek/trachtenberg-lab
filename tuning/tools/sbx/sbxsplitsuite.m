function sbxsplitsuite(fout)

% splits fout into the individual sessions

load('-mat',fout);  % load matlab file for merged
% S = load('-mat',[fout '_suite.signals']); % load signals
S = load('-mat',[fout '.signals']); % load signals

nf = length(info.source_max_idx);   % # of files
filepath = strcat(fileparts(fout), filesep);

k = 0;
for i = 1:nf
    idx  = (k+1) : (k+info.source_max_idx(i));
    sig  = S.sig(idx,:);
    np   = S.np(idx,:);
    spks = S.spks(idx,:);
    save([filepath, info.source{i} '_suite.signals'],'sig','np','spks');
    k = idx(end);
end

