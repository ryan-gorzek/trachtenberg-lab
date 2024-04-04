
addpath(genpath(fullfile("E:/Imaging/code/tools/")));

%% merge sessions

% move .sbx and .mat files out into current directory

sbxmergesessions('pv_vip_06_444_444', {'pv_vip_06_410_000', ...
                                       'pv_vip_06_420_000'});

% then run suite2p

%% convert suite2p segmentation to sbx

sbxsuite2sbx('pv_vip_06_444_444/suite2p/plane0/Fall.mat', 'E:/Imaging/pv_vip_06/pv_vip_06_444_444/pv_vip_06_444_444');

%% generate .signals file

sbxf2spks('E:/Imaging/pv_vip_06/pv_vip_06_444_444/pv_vip_06_444_444');

%% generate matfile for splitting .signals

sbxsplitsuite('E:/Imaging/pv_vip_06/pv_vip_06_444_444/pv_vip_06_444_444');
