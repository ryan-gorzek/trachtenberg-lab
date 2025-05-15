
%% Mouse 3D

% Load data from CSV
data = readtable('E:/Transcriptomics_V1/Integration/PCs/IT/mouse_pc_data_shared_subsample_within_with_subclass.csv');

% Extract the first three PCs
pc1 = data.pca_1;  % Assuming first column is PC1
pc2 = data.pca_2;  % Assuming second column is PC2
pc3 = data.pca_3 * -1;  % Assuming third column is PC3

% Convert subclass to categorical and assign unique numerical values
[subclass_groups, ~, subclass_idx] = unique(data.subclass, 'stable');

% Define custom hex colors
hex_colors = {'#FFB3B3'; '#FF7F50'; '#FFA07A'; '#FF6347'}; % Example hex colors

% Convert hex colors to RGB (normalize to [0,1] range)
cmap = cell2mat(cellfun(@(x) sscanf(x(2:end),'%2x%2x%2x')'/255, hex_colors, 'UniformOutput', false));

% Create 3D scatter plot
figure;
set(gcf,'Renderer','Painters')
hold on;
for i = 1:length(subclass_groups)
    idx = (subclass_idx == i);
    scatter3(pc1(idx), pc2(idx), pc3(idx), 15, cmap(i, :), 'filled', 'DisplayName', subclass_groups{i});
end

% Load the tetrahedron coordinates from MAT file
load('Mouse_IT_Shared_Within_30PCs.mat', 'arcOrig');

% Plot the tetrahedron edges
% Assuming 'arc' contains a 4x3 matrix where each row is a vertex [x, y, z]
tetra_edges = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; % Indices defining the edges of a tetrahedron
for k = 1:size(tetra_edges, 1)
    plot3(arcOrig(tetra_edges(k, :), 1), arcOrig(tetra_edges(k, :), 2), arcOrig(tetra_edges(k, :), 3) * -1, '-k', 'LineWidth', 1.5);
end

axis equal;
axis square;
view(156, 18);
hold off;

xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title('3D Scatter Plot of Principal Components');
grid on;

view(156, 18);
% view(-12.7, 18);

%% Opossum 3D

% Load data from CSV
data = readtable('E:/Transcriptomics_V1/Integration/PCs/IT/opossum_pc_data_shared_subsample_within_with_subclass.csv');

% Extract the first three PCs
pc1 = data.pca_1 * -1;  % Assuming first column is PC1
pc2 = data.pca_3;  % Assuming second column is PC2
pc3 = data.pca_2;  % Assuming third column is PC3

% Convert subclass to categorical and assign unique numerical values
[subclass_groups, ~, subclass_idx] = unique(data.subclass, 'stable');

% Define custom hex colors
hex_colors = {'#FFB3B3'; '#FF7F50'; '#FF6347'; '#FFA07A'}; % Example hex colors

% Convert hex colors to RGB (normalize to [0,1] range)
cmap = cell2mat(cellfun(@(x) sscanf(x(2:end),'%2x%2x%2x')'/255, hex_colors, 'UniformOutput', false));

% Create 3D scatter plot
figure;
set(gcf,'Renderer','Painters')
hold on;
for i = 1:length(subclass_groups)
    idx = (subclass_idx == i);
    scatter3(pc1(idx), pc2(idx), pc3(idx), 15, cmap(i, :), 'filled', 'DisplayName', subclass_groups{i});
end

% Load the tetrahedron coordinates from MAT file
load('Opossum_IT_Shared_Within_30PCs.mat', 'arcOrig');

% Plot the tetrahedron edges
% Assuming 'arc' contains a 4x3 matrix where each row is a vertex [x, y, z]
tetra_edges = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; % Indices defining the edges of a tetrahedron
for k = 1:size(tetra_edges, 1)
    plot3(arcOrig(tetra_edges(k, :), 1) * -1, arcOrig(tetra_edges(k, :), 3), arcOrig(tetra_edges(k, :), 2), '-k', 'LineWidth', 1.5);
end

axis equal;
axis square;
view(275, 25);
hold off;

xlabel('PC1');
ylabel('PC3');
zlabel('PC2');
title('3D Scatter Plot of Principal Components');
grid on;
% legend('show', 'Location', 'best');
