
addpath(genpath(fullfile("E:/code/tools/")));

%% set parameters

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
             sprintf('E:/%s/%s/%s', animal, merge_id, merge_id));

%% generate .signals file

sbxf2spks(sprintf('E:/%s/%s/%s', animal, merge_id, merge_id));

%% generate matfile for splitting .signals

sbxsplitsuite(sprintf('E:/%s/%s/%s', animal, merge_id, merge_id));

%% archive

% % merge separate sessions for battery3 and battery4
% 
% ftypes = {'.mat', '.sbx', '_p.mat', '_eye.mat', '_quadrature.mat'};
% 
% if any(order == 5)
%     segs = strcat(animal, '_', n, {'50_001', '50_002', '50_003'});
%     newfile = sprintf('%s_%s50_000', animal, n);
%     for s = segs
%         for f = ftypes
%             movefile([s{1}, filesep, s{1}, f{1}], pwd);
%         end
%     end
%     % merge scanbox and metadata
%     sbxmergesessions(newfile, segs);
%     % merge stimulus data
%     mergeparams(newfile, segs);
%     % merge eye data
%     mergeeye(newfile, segs);
%     % merge quad data
%     mergequad(newfile, segs);
%     % move all the files back and make new folder for merged recording
%     for s = segs
%         for f = ftypes
%             movefile([s{1}, f{1}], [s{1}, filesep]);
%         end
%     end
%     mkdir(newfile);
%     for f = ftypes
%         movefile([newfile, f{1}], [newfile, filesep]);
%     end
% end
% 
% function mergeparams(newfile, mergefiles)
%     param_all = [];
%     perms_all = [];
%     for f = mergefiles
%         load([f{1}, '_p.mat'], 'param', 'perms');
%         param_all = param;
%         perms_all = horzcat(perms_all, perms);
%     end
%     param = param_all;
%     perms = perms_all;
%     save([newfile, '_p.mat'], 'param', 'perms');
% end
% 
% function mergeeye(newfile, mergefiles)
%     abstime_all = [];
%     data_all = [];
%     time_all = [];
%     for f = mergefiles
%         load([f{1}, '_eye.mat'], 'abstime', 'data', 'time');
%         abstime_all = vertcat(abstime_all, abstime);
%         data_all = cat(4, data_all, data);
%         time_all = vertcat(time_all, time);
%     end
%     abstime = abstime_all;
%     data = data_all;
%     time = time_all;
%     save([newfile, '_eye.mat'], 'abstime', 'data', 'time');
% end
% 
% function mergequad(newfile, mergefiles)
%     quad_data_all = [];
%     for f = mergefiles
%         load([f{1}, '_quadrature.mat'], 'quad_data');
%         quad_data_all = horzcat(quad_data_all, quad_data);
%     end
%     quad_data = quad_data_all;
%     save([newfile, '_quadrature.mat'], 'quad_data');
% end
