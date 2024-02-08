function sbxsplitsessions(fout)

% splits fout into the individual sessions

load('-mat',fout);  % load matlab file for merged
S = load('-mat',[fout '_rigid.signals']); % load signals

nf = length(info.source_max_idx);   % # of files

k = 1;
for i = 1:nf
    idx  = k : (k+info.source_max_idx(i));
    sig  = S.sig(idx,:);
    np   = S.np(idx,:);
    spks = S.spks(idx,:);
    save([info.source{i} '_rigid.signals'],'sig','np','spks');
    k = idx(end)+1;
end

