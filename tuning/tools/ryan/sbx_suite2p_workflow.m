
addpath(genpath(fullfile("E:/Imaging/code/tools/")));

%% merge sessions

% move .sbx and .mat files out into current directory

sbxmergesessions('pv_vip_01_444_444', {'pv_vip_01_410_000', ...
                                       'pv_vip_01_420_000', ...
                                       'pv_vip_01_430_000', ...
                                       'pv_vip_01_440_000', ...
                                       'pv_vip_01_450_000', ...
                                       'pv_vip_01_460_000'});

% then run suite2p

%% convert suite2p segmentation to sbx

sbxsuite2sbx('pv_vip_01_444_444/suite2p/plane0/Fall.mat', 'E:/Imaging/pv_vip_01/pv_vip_01_444_444/pv_vip_01_444_444');

%% generate .signals file

sbxf2spks('E:/Imaging/pv_vip_01/pv_vip_01_444_444/pv_vip_01_444_444');

%% generate matfile for splitting .signals

sbxsplitsuite('E:/Imaging/pv_vip_01/pv_vip_01_444_444/pv_vip_01_444_444');
