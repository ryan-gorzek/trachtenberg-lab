
%% Mouse (Shared & Subsampled Space)

addpath(genpath(fullfile("C:/Users/TLab/Documents/Ryan/particode/")));

% Load data from CSV
input = readtable('E:/Transcriptomics_V1/Integration/PCs/IT/mouse_pc_data_shared_subsample_with_subclass.csv');

% Extract the first three PCs
data = [];
for pc = 1:30
    data = [data, input.(sprintf('pca_%i', pc))];  % Assuming first column is PC1
end

labels = input.subclass;
ident = {'subclass'};

[arc, arcOrig, pc, errs, pval, coefs] = ParTI(data, 1, 5, labels, ...
    ident, 0, [], [], [], 0.2, 'Mouse_IT_Shared_30PCs');

%% Opossum (Shared & Subsampled Space)

addpath(genpath(fullfile("C:/Users/TLab/Documents/Ryan/particode/")));

% Load data from CSV
input = readtable('E:/Transcriptomics_V1/Integration/PCs/IT/opossum_pc_data_shared_subsample_within_with_subclass.csv');

% Extract the first three PCs
data = [];
for pc = 1:30
    data = [data, input.(sprintf('pca_%i', pc))];  % Assuming first column is PC1
end

labels = input.subclass;
ident = {'subclass'};

[arc, arcOrig, pc, errs, pval, coefs] = ParTI(data, 1, 5, labels, ...
    ident, 0, [], [], [], 0.2, 'Opossum_IT_Shared_30PCs');
