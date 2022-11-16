
%% Preprocessing

clearvars; close all;
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

clearvars; close all;
file_path = "/Users/ryan.gorzek/Library/CloudStorage/GoogleDrive-ryan.gorzek@gmail.com/Shared drives/Astrocytes - Protocadherins/";
cells_tbl = load(strcat(file_path, "randorisf_cells_tbl.mat")).cells_tbl;

subj_ids = ["john01", "john02", "john03"]; % 
genotypes = ["KO", "KO", "Het"];
plane_ids = ["000", "001", "002"]; % 
for subj = subj_ids
    % First, plot tuning kernels for each eye from all planes in each animal.
    subj_idx = cells_tbl.subject_id == subj;
    subj_tbl = cells_tbl(subj_idx, :);
    title = sprintf("%s (%s) %i cells", subj, genotypes(subj_ids == subj), size(subj_tbl, 1));
    save_path = strcat(file_path, "plots", filesep, subj, "_all_kernels.jpeg");
    kernPlot(subj_tbl, num_rows=15, title=title, save_path=save_path);
    %
    for plane = plane_ids
        % Then, plot tuning kernels for each eye from each plane.
        plane_idx = subj_idx & cells_tbl.plane_id == plane;
        plane_tbl = cells_tbl(plane_idx, :);
        title = sprintf("%s_%s (%s) %i cells", subj, plane, genotypes(subj_ids == subj), size(plane_tbl, 1));
        save_path = strcat(file_path, "plots", filesep, subj, "_", plane, "_kernels.jpeg");
        kernPlot(plane_tbl, title=title, save_path=save_path);
        %
    end
end

% Define plotting functions.

function kernPlot(cells_table, NVAs)
% KERNPLOT Plot orisf tuning kernels from a table of sbxorisf_new data.
arguments
    cells_table
    NVAs.num_rows (1,1) double = 10
    NVAs.num_cols (1,1) double = NaN
    NVAs.title string = ""
    NVAs.save_path string = ""
end

num_cells = size(cells_table, 1);
if isnan(NVAs.num_cols)
    NVAs.num_cols = ceil(num_cells/NVAs.num_rows);
end
% Create a 5x10 tile area for each cell, spacing columns by 1 tile.
f = figure; hold on; set(gcf, "Color","w");
t = tiledlayout(NVAs.num_rows*5, ...
                NVAs.num_cols*10, ...
                "TileSpacing","tight", ...
                "Padding","compact");
for cell_num = 1:num_cells
    % Concatenate ipsi and contra side-by-side, separated by a NaN column.
    kern_contra = cells_table{cell_num, "kern_contra"}{:};
    kern_ipsi = cells_table{cell_num, "kern_ipsi"}{:};
    kerns = horzcat(kern_ipsi, nan(size(kern_contra,1), 1), kern_contra);
    % Plot.
    nexttile([5,10]);
    imagesc(kerns, "AlphaData",~isnan(kerns));
    axis tight; axis off; pbaspect([2.04,1,1]);
    if NVAs.num_rows == 10 && ...
       cell_num <= NVAs.num_cols
        title("  I              C ");
    elseif cell_num <= NVAs.num_cols
        title(" I       C");
    end
end
% Dynamically resize the figure such that the kernel plots appear similar
% across animals/planes.
if NVAs.num_rows == 10
    f.Position = [-231, 1130, 125*NVAs.num_cols, 70*NVAs.num_rows];
else
    f.Position = [-231, 1130, 78*NVAs.num_cols, 50*NVAs.num_rows];
end
% Set title and axis labels.
title(t, NVAs.title, ...
         "Interpreter", "none", ...
         "FontWeight", "bold", ...
         "FontSize", 15);
annotation("textbox", ...
           "String","Orientation", ...
           "FontWeight","bold", ...
           "FontSize",15, ...
           "EdgeColor","none", ...
           "HorizontalAlignment","center", ...
           "VerticalAlignment","middle", ...
           "Rotation",90, ...
           "Position",[0.02, 0.5, 0, 0]);
annotation("textbox", ...
           "String","Spatial Frequency", ...
           "FontWeight","bold", ...
           "FontSize",15, ...
           "EdgeColor","none", ...
           "HorizontalAlignment","center", ...
           "VerticalAlignment","middle", ...
           "Rotation",0, ...
           "Position",[0.4, 0.025, 0.2, 0]);
% Set one y-axis and one x-axis per plot for reference.

% Save the figure to file.
imwrite(getframe(f).cdata, ...
        NVAs.save_path, ...
        "jpeg", ...
        "Mode", "lossless", ...
        "Quality", 100);

end
