
%% Mouse (Original PC Space)

addpath(genpath(fullfile("C:/Users/TLab/Documents/Ryan/particode/")));

% Load data from CSV
input = readtable('E:/Transcriptomics_V1/Integration/PCs/IT/mouse_pc_data_subsample_with_subclass.csv');

% Extract the first three PCs
data = [];
for pc = 1:10
    data = [data, input.(sprintf('pca_%i', pc))];  % Assuming first column is PC1
end

labels = input.subclass;
ident = {'subclass'};

[arc, arcOrig, pc, errs, pval, coefs] = ParTI(data, 1, 5, labels, ...
    ident, 0, [], [], [], 0.2, 'Mouse_IT_3PCs');

%% Mouse (GE Space)

addpath(genpath(fullfile("C:/Users/TLab/Documents/Ryan/particode/")));

% Load data from CSV
input = readtable('E:/Transcriptomics_V1/Integration/PCs/IT/mouse_ge_data_subsample_with_subclass.csv');

data = table2array(input(:, 1:end-1));

labels = input.subclass;
ident = {'subclass'};

[arc, arcOrig, pc, errs, pval, coefs] = ParTI(data, 1, 5, labels, ...
    ident, 0, [], [], [], 0.2, 'Mouse_IT_GE');

%% Opossum (Original Space)

addpath(genpath(fullfile("C:/Users/TLab/Documents/Ryan/particode/")));

% Load data from CSV
input = readtable('E:/opossum_pc_data_subsample_with_subclass.csv');

% Extract the first three PCs
data = [];
for pc = 1:30
    data = [data, input.(sprintf('pca_%i', pc))];  % Assuming first column is PC1
end

labels = input.subclass;
ident = {'subclass'};

[arc, arcOrig, pc, errs, pval, coefs] = ParTI(data, 1, 5, labels, ...
    ident, 0, [], [], [], 0.2, 'Opossum_IT_3PCs');

%% Opossum (GE Space)

addpath(genpath(fullfile("C:/Users/TLab/Documents/Ryan/particode/")));

% Load data from CSV
input = readtable('E:/Transcriptomics_V1/Integration/PCs/IT/opossum_ge_data_subsample_with_subclass.csv');

data = table2array(input(:, 1:end-1));

labels = input.subclass;
ident = {'subclass'};

[arc, arcOrig, pc, errs, pval, coefs] = ParTI(data, 1, 5, labels, ...
    ident, 0, [], [], [], 0.2, 'Opossum_IT_GE');

%% Opossum (Original Space without Projection)

addpath(genpath(fullfile("C:/Users/TLab/Documents/Ryan/particode/")));

% Load data from CSV
input = readtable('E:/opossum_pc_data_subsample_noproject_with_subclass.csv');

% Extract the first three PCs
data = [];
for pc = 1:30
    data = [data, input.(sprintf('PC_%i', pc))];  % Assuming first column is PC1
end

labels = input.subclass;
ident = {'subclass'};

[arc, arcOrig, pc, errs, pval, coefs] = ParTI(data, 1, 5, labels, ...
    ident, 0, [], [], [], 0.2, 'Opossum_IT_NoProject_3PCs');

%% Opossum (Shared & Subsampled Space)

addpath(genpath(fullfile("C:/Users/TLab/Documents/Ryan/particode/")));

% Load data from CSV
input = readtable('E:/opossum_pc_data_shared_subsample_with_subclass.csv');

% Extract the first three PCs
data = [];
for pc = 1:30
    data = [data, input.(sprintf('pca_%i', pc))];  % Assuming first column is PC1
end

labels = input.subclass;
ident = {'subclass'};

[arc, arcOrig, pc, errs, pval, coefs] = ParTI(data, 1, 5, labels, ...
    ident, 0, [], [], [], 0.2, 'Opossum_IT_Shared_30PCs');
