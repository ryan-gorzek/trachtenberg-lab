
%% Shared (Subsampled, Project onto Mouse)

% Load data from CSV
data = readtable('E:/Transcriptomics_V1/Integration/PCs/IT/shared_mouse_pc_data_downsample_subsample_with_subclass.csv');
% data = data(ismember(data.species, 'Mouse'), :);
data(ismember(data.species, 'Mouse'), "subclass") = repmat({'Mouse'}, [nnz(ismember(data.species, 'Mouse')), 1]);

% Extract the first three PCs
pc1 = data.pca_1;  % Assuming first column is PC1
pc2 = data.pca_2;  % Assuming second column is PC2
pc3 = data.pca_3 * -1;  % Assuming third column is PC3

% Convert subclass to categorical and assign unique numerical values
[subclass_groups, ~, subclass_idx] = unique(data.subclass, 'stable');

% Define a colormap based on the number of unique subclasses
hex_colors = {'#AAAAAA'; '#FFB3B3'; '#FFA07A'; '#FF6347'; '#FF7F50'}; % Example hex colors

% Convert hex colors to RGB (normalize to [0,1] range)
cmap = cell2mat(cellfun(@(x) sscanf(x(2:end),'%2x%2x%2x')'/255, hex_colors, 'UniformOutput', false));

% Create 3D scatter plot
figure;
set(gcf,'Renderer','Painters')
hold on;
for i = length(subclass_groups):-1:1
    idx = (subclass_idx == i);
    scatter3(pc1(idx), pc2(idx), pc3(idx), 5, cmap(i, :), 'filled', 'DisplayName', subclass_groups{i});
end

% % Load the tetrahedron coordinates from MAT file
% load('Mouse_IT_3PCs_Workspace.mat', 'arcOrig');

% % Plot the tetrahedron edges
% % Assuming 'arc' contains a 4x3 matrix where each row is a vertex [x, y, z]
% tetra_edges = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; % Indices defining the edges of a tetrahedron
% for k = 1:size(tetra_edges, 1)
%     plot3(arcOrig(tetra_edges(k, :), 1) * 0.75, arcOrig(tetra_edges(k, :), 2) * 0.75, arcOrig(tetra_edges(k, :), 3) * 0.75, '-k', 'LineWidth', 1.5);
% end

axis equal;
axis square;
view(-25, 22);
hold off;

xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title('3D Scatter Plot of Principal Components');
grid on;

%% Shared (Subsampled, Project onto Opossum)

% Load data from CSV
data = readtable('E:/Transcriptomics_V1/Integration/PCs/IT/shared_opossum_pc_data_downsample_subsample_with_subclass.csv');
% data = data(ismember(data.species, 'Mouse'), :);
data(ismember(data.species, 'Opossum'), "subclass") = repmat({'Opossum'}, [nnz(ismember(data.species, 'Opossum')), 1]);

% Extract the first three PCs
pc1 = data.pca_1;  % Assuming first column is PC1
pc2 = data.pca_2;  % Assuming second column is PC2
pc3 = data.pca_3 * -1;  % Assuming third column is PC3

% Convert subclass to categorical and assign unique numerical values
[subclass_groups, ~, subclass_idx] = unique(data.subclass, 'stable');

% Define a colormap based on the number of unique subclasses
hex_colors = {'#FFB3B3'; '#FFA07A'; '#FF6347'; '#FF7F50'; '#AAAAAA'}; % Example hex colors

% Convert hex colors to RGB (normalize to [0,1] range)
cmap = cell2mat(cellfun(@(x) sscanf(x(2:end),'%2x%2x%2x')'/255, hex_colors, 'UniformOutput', false));

% Create 3D scatter plot
figure;
set(gcf,'Renderer','Painters')
hold on;
for i = length(subclass_groups):-1:1
    idx = (subclass_idx == i);
    scatter3(pc1(idx), pc2(idx), pc3(idx), 5, cmap(i, :), 'filled', 'DisplayName', subclass_groups{i});
end

% % Load the tetrahedron coordinates from MAT file
% load('Mouse_IT_3PCs_Workspace.mat', 'arcOrig');

% % Plot the tetrahedron edges
% % Assuming 'arc' contains a 4x3 matrix where each row is a vertex [x, y, z]
% tetra_edges = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; % Indices defining the edges of a tetrahedron
% for k = 1:size(tetra_edges, 1)
%     plot3(arcOrig(tetra_edges(k, :), 1) * 0.75, arcOrig(tetra_edges(k, :), 2) * 0.75, arcOrig(tetra_edges(k, :), 3) * 0.75, '-k', 'LineWidth', 1.5);
% end

axis equal;
axis square;
view(-25, 22);
hold off;

xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title('3D Scatter Plot of Principal Components');
grid on;
