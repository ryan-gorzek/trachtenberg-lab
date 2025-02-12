
%% Heatmaps (Euclidean Distance)

% Load data from CSV for both species
mouse_data = readtable('E:/mouse_pc_data_subsample_with_subclass.csv');
opossum_data = readtable('E:/opossum_pc_data_subsample_with_subclass.csv');

% Extract the first three PCs for both species
mouse_pc_matrix = [mouse_data.pca_1, mouse_data.pca_2, mouse_data.pca_3];
opossum_pc_matrix = [opossum_data.pca_1, opossum_data.pca_2, opossum_data.pca_3];

% Convert subclass to categorical and assign numerical labels
mouse_subclasses = unique(mouse_data.subclass);
opossum_subclasses = unique(opossum_data.subclass);

% Initialize matrices to store Euclidean distances for subclass pairs
mouse_distance_matrix = zeros(length(mouse_subclasses));
opossum_distance_matrix = zeros(length(opossum_subclasses));

% Compute Euclidean distances for each subclass pair in mouse
for i = 1:length(mouse_subclasses)
    for j = 1:length(mouse_subclasses)
        if i ~= j  % Avoid diagonal computation
            % Select cells belonging to the two subclasses
            idx_i = strcmp(mouse_data.subclass, mouse_subclasses{i});
            idx_j = strcmp(mouse_data.subclass, mouse_subclasses{j});
            
            % Compute centroid of each subclass in PC space
            centroid_i = mean(mouse_pc_matrix(idx_i, :), 1);
            centroid_j = mean(mouse_pc_matrix(idx_j, :), 1);
            
            % Compute Euclidean distance
            mouse_distance_matrix(i, j) = norm(centroid_i - centroid_j);
        end
    end
end

% Compute Euclidean distances for each subclass pair in opossum
for i = 1:length(opossum_subclasses)
    for j = 1:length(opossum_subclasses)
        if i ~= j  % Avoid diagonal computation
            % Select cells belonging to the two subclasses
            idx_i = strcmp(opossum_data.subclass, opossum_subclasses{i});
            idx_j = strcmp(opossum_data.subclass, opossum_subclasses{j});
            
            % Compute centroid of each subclass in PC space
            centroid_i = mean(opossum_pc_matrix(idx_i, :), 1);
            centroid_j = mean(opossum_pc_matrix(idx_j, :), 1);
            
            % Compute Euclidean distance
            opossum_distance_matrix(i, j) = norm(centroid_i - centroid_j);
        end
    end
end

% Set diagonal to zero to appear white
mouse_distance_matrix(eye(size(mouse_distance_matrix)) == 1) = 0;
opossum_distance_matrix(eye(size(opossum_distance_matrix)) == 1) = 0;

% Set upper right triangle to zero for mouse and lower left for opossum
mouse_distance_matrix(triu(true(size(mouse_distance_matrix)), 1)) = 0;
opossum_distance_matrix(tril(true(size(opossum_distance_matrix)), -1)) = 0;

% Define color scale limits dynamically based on the range of distances
clims = [min([mouse_distance_matrix(:); opossum_distance_matrix(:)]), ...
         max([mouse_distance_matrix(:); opossum_distance_matrix(:)])];

% Define colormap (white to red to blue)
custom_colormap = [linspace(1,1,256)', linspace(1,0,256)', linspace(1,0,256)';  % White to red
                   linspace(1,0,256)', zeros(256,1), linspace(0,1,256)'];       % Red to blue

% Plot heatmap for mouse
figure;
h_mouse = heatmap(mouse_subclasses, mouse_subclasses, round(mouse_distance_matrix,2), ...
    'Colormap', custom_colormap, 'ColorLimits', clims);
title('Mouse IT Subclass Euclidean Distance Heatmap');
xlabel('Subclass');
ylabel('Subclass');
h_mouse.CellLabelFormat = '%.2f';  % Round values to 2 decimal places
h_mouse.GridVisible = 'off';  % Remove gridlines for a cleaner look

% Plot heatmap for opossum
figure;
h_opossum = heatmap(opossum_subclasses, opossum_subclasses, round(opossum_distance_matrix,2), ...
    'Colormap', custom_colormap, 'ColorLimits', clims);
title('Opossum IT Subclass Euclidean Distance Heatmap');
xlabel('Subclass');
ylabel('Subclass');
h_opossum.CellLabelFormat = '%.2f';  % Round values to 2 decimal places
h_opossum.GridVisible = 'off';  % Remove gridlines for a cleaner look
