
addpath(genpath(fullfile("E:/code/tools/")));

%% set parameters

drive = "E:";
animal = 'pv_vip_06';
n = '7';
order = [1, 2, 5, 6, 3, 4];

merge_id = sprintf('%s_%s%s%s_%s%s%s', animal, n, n, n, n, n, n);

%% merge all stimuli

ftypes = {'.mat', '.sbx'};
sessions = cell(1, numel(order));
for s = 1:numel(order)
    sessions{s} = sprintf('%s_%s%s0_000', animal, n, num2str(order(s)));
    for f = ftypes
        movefile([sessions{s}, filesep, sessions{s}, f{1}], pwd);
    end
end
sbxmergesessions(merge_id, sessions);
mkdir(merge_id);
sessions = [sessions, {merge_id}];
for s = sessions
    for f = ftypes
        movefile([s{1}, f{1}], [s{1}, filesep]);
    end
end

% then run suite2p

%% convert suite2p segmentation to sbx

sbxsuite2sbx(sprintf('%s/suite2p/plane0/Fall.mat', merge_id), ...
             sprintf('%s/%s/%s/%s', drive, animal, merge_id, merge_id));

%% generate .signals file

sbxf2spks(sprintf('%s/%s/%s/%s', drive, animal, merge_id, merge_id));

%% split .signals into original sessions

sbxsplitsuite(sprintf('%s/%s/%s/%s', drive, animal, merge_id, merge_id));
for s = 1:numel(order)
    session = sprintf('%s_%s%s0_000', animal, n, num2str(order(s)));
    movefile([merge_id, filsep, session, '_suite.signals'], session); %%%% see if this works
end
