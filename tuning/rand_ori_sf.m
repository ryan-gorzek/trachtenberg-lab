
%% Preprocessing

file_path = "/Users/ryan.gorzek/Library/CloudStorage/GoogleDrive-ryan.gorzek@gmail.com/Shared drives/Astrocytes - Protocadherins/";
addpath(genpath(fullfile(file_path)));

subj_ids = ["john01", "john02", "john03"];
genotypes = ["KO", "KO", "Het"];
plane_ids = ["000", "001", "002"];
eye_ids = ["000", "001"];
eye_names = ["contra", "ipsi"];
% Preallocate table for storing cell stats across subjects and planes.
cells_tbl = [];
for subj = subj_ids
    for plane = plane_ids
        exp_id = char(strcat(subj, "_", plane));
        fprintf("Processing %s...\n", exp_id);
        exp_path = char(strcat(file_path, filesep, subj, filesep, exp_id, filesep));
        plane_path = strcat(exp_path, "suite2p", filesep, "plane0", filesep);
        % Generate .align, .segment, and .signals1 files. 
        sbxsuite2sbx(strcat(plane_path, "Fall"), strcat(exp_path, exp_id));
        % Generate .signals file.
        sbxf2spks(strcat(exp_path, exp_id));
        % Split .signals file for the two eyes.
        sbxsplitsuite(strcat(exp_path, exp_id));
        % Process contra and ipsi eyes.
        eye_tbls = cell(1,2);
        for eyeball = eye_ids
            % Calculate tuning results and generate .orisf file for each eye.
            eye_id = char(strcat(exp_path, exp_id, "_", eyeball, "_suite"));
            sbxorisf_new(eye_id);
            % Load .orisf file and create table from structure array.
            stat_struct = load(strcat(eye_id, ".orisf"), "-mat", "stat").stat;
            eye_idx = eye_ids == eyeball;
            eye_tbls{eye_idx} = struct2table(stat_struct);
            % Add eye labels to stat table variable names.
            eye_tbls{eye_idx}.Properties.VariableNames = ...
                strcat(eye_tbls{eye_idx}.Properties.VariableNames, "_", eye_names{eye_idx});
        end
        % Add table variable for subject and plane IDs.
        plane_tbl = [eye_tbls{1}, eye_tbls{2}];
        num_cells = size(plane_tbl, 1);
        plane_id = repmat(plane, [num_cells, 1]);
        plane_tbl = addvars(plane_tbl, plane_id, 'Before', "k_contra");
        genotype = repmat(genotypes(subj_ids == subj), [num_cells, 1]);
        plane_tbl = addvars(plane_tbl, genotype, 'Before', "plane_id");
        subject_id = repmat(subj, [num_cells, 1]);
        plane_tbl = addvars(plane_tbl, subject_id, 'Before', "genotype");
        % Add current data to cumulative table.
        cells_tbl = vertcat(cells_tbl, plane_tbl);
    end
end

% Save table containing data for all subjects and planes.
save(strcat(file_path, "randorisf_cells_tbl.mat"), "cells_tbl");

%% Plotting

% Plot kernels for each eye.
