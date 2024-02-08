function sbxmergesessions(fout,fn)

% merges different experiments into a single file
% fn is a cell array of file names
% fout is the new, merged data

global info

nf = length(fn);

% first store the infos

cmd = ['copy /b '];
max_idx = [];
for i = 1:length(fn)
    z = sbxread(fn{i},0,1);
    max_idx = [max_idx info.max_idx];
    cmd = [cmd fn{i} '.sbx+'];
end

info.source = fn;
info.source_max_idx = max_idx;

if ~exist([fout '.sbx'],'file')
    cmd = [cmd(1:end-1) ' ' fout '.sbx'];
    disp('Please wait...');
    system(cmd);
end

if ~exist([fout '.mat'],'file')
    save([fout '.mat'],'info');
end


