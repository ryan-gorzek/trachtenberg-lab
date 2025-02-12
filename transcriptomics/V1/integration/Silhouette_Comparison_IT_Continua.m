
%% Heatmaps (Silhouette)

% Load data from CSV for both species
mouse_data = readtable('E:/mouse_pc_data_subsample_with_subclass.csv');
opossum_data = readtable('E:/opossum_pc_data_subsample_with_subclass.csv');

% Extract the first three PCs for both species
mouse_pc_matrix = [mouse_data.pca_1, mouse_data.pca_2, mouse_data.pca_3];
opossum_pc_matrix = [opossum_data.pca_1, opossum_data.pca_2, opossum_data.pca_3];

% Convert subclass to categorical and assign numerical labels
mouse_subclasses = unique(mouse_data.subclass);
opossum_subclasses = unique(opossum_data.subclass);

% Initialize matrices to store silhouette scores for subclass pairs
mouse_silhouette_matrix = zeros(length(mouse_subclasses));
opossum_silhouette_matrix = zeros(length(opossum_subclasses));

% Compute silhouette scores for each subclass pair in mouse
for i = 1:length(mouse_subclasses)
    for j = 1:length(mouse_subclasses)
        if i ~= j  % Avoid diagonal computation
            % Select cells belonging to the two subclasses
            idx = strcmp(mouse_data.subclass, mouse_subclasses{i}) | strcmp(mouse_data.subclass, mouse_subclasses{j});
            subclass_labels = mouse_data.subclass(idx);
            pc_subset = mouse_pc_matrix(idx, :);
            
            % Convert labels to numerical indices
            numeric_labels = grp2idx(categorical(subclass_labels));
            
            % Compute silhouette score
            silhouette_scores = silhouette(pc_subset, numeric_labels);
            mouse_silhouette_matrix(i, j) = median(silhouette_scores);
        end
    end
end

% Compute silhouette scores for each subclass pair in opossum
for i = 1:length(opossum_subclasses)
    for j = 1:length(opossum_subclasses)
        if i ~= j  % Avoid diagonal computation
            % Select cells belonging to the two subclasses
            idx = strcmp(opossum_data.subclass, opossum_subclasses{i}) | strcmp(opossum_data.subclass, opossum_subclasses{j});
            subclass_labels = opossum_data.subclass(idx);
            pc_subset = opossum_pc_matrix(idx, :);
            
            % Convert labels to numerical indices
            numeric_labels = grp2idx(categorical(subclass_labels));
            
            % Compute silhouette score
            silhouette_scores = silhouette(pc_subset, numeric_labels);
            opossum_silhouette_matrix(i, j) = median(silhouette_scores);
        end
    end
end

% Set diagonal to zero to appear white
mouse_silhouette_matrix(eye(size(mouse_silhouette_matrix)) == 1) = 0;
opossum_silhouette_matrix(eye(size(opossum_silhouette_matrix)) == 1) = 0;

% Set upper right triangle to zero for mouse and lower left for opossum
mouse_silhouette_matrix(triu(true(size(mouse_silhouette_matrix)), 1)) = 0;
opossum_silhouette_matrix(tril(true(size(opossum_silhouette_matrix)), -1)) = 0;

% Set color scale limits from 0.5 to 1
clims = [0.5, 1];

% Define white-red-blue colormap
custom_colormap = [linspace(1,1,256)', linspace(1,0,256)', linspace(1,0,256)';  % White to red
                   linspace(1,0,256)', zeros(256,1), linspace(0,1,256)'];       % Red to blue

% Plot heatmap for mouse
figure;
h_mouse = heatmap(mouse_subclasses, mouse_subclasses, round(mouse_silhouette_matrix,2), ...
    'Colormap', custom_colormap, 'ColorLimits', clims);
title('Mouse IT Subclass Silhouette Score Heatmap');
xlabel('Subclass');
ylabel('Subclass');
h_mouse.CellLabelFormat = '%.2f';  % Round values to 2 decimal places
h_mouse.GridVisible = 'off';  % Remove gridlines for a cleaner look

% Plot heatmap for opossum
figure;
h_opossum = heatmap(opossum_subclasses, opossum_subclasses, round(opossum_silhouette_matrix,2), ...
    'Colormap', custom_colormap, 'ColorLimits', clims);
title('Opossum IT Subclass Silhouette Score Heatmap');
xlabel('Subclass');
ylabel('Subclass');
h_opossum.CellLabelFormat = '%.2f';  % Round values to 2 decimal places
h_opossum.GridVisible = 'off';  % Remove gridlines for a cleaner look

%% Heatmaps (Cross-Species Silhouette Comparison)

% Load data from CSV for both species projected into the same PC space
mouse_data = readtable('E:/mouse_pc_data_shared_subsample_with_subclass.csv');
opossum_data = readtable('E:/opossum_pc_data_shared_subsample_with_subclass.csv');

% Extract the first three PCs for both species
mouse_pc_matrix = [mouse_data.pca_1, mouse_data.pca_2, mouse_data.pca_3];
opossum_pc_matrix = [opossum_data.pca_1, opossum_data.pca_2, opossum_data.pca_3];

% Convert subclass to categorical and assign numerical labels
mouse_subclasses = unique(mouse_data.subclass);
opossum_subclasses = unique(opossum_data.subclass);

% Combine data for cross-species comparison
combined_pc_matrix = [mouse_pc_matrix; opossum_pc_matrix];
combined_subclass_labels = [mouse_data.subclass; opossum_data.subclass];

% Get unique subclass names across both species
combined_subclasses = unique(combined_subclass_labels);

% Initialize matrix to store silhouette scores for cross-species subclass pairs
cross_species_silhouette_matrix = zeros(length(mouse_subclasses), length(opossum_subclasses));

% Compute silhouette scores for each subclass pair across species
for i = 1:length(mouse_subclasses)
    for j = 1:length(opossum_subclasses)
        % Select cells belonging to the mouse and opossum subclasses
        idx = strcmp(combined_subclass_labels, mouse_subclasses{i}) | strcmp(combined_subclass_labels, opossum_subclasses{j});
        subclass_labels_subset = combined_subclass_labels(idx);
        pc_subset = combined_pc_matrix(idx, :);
        
        % Convert subclass labels to numerical indices
        numeric_labels = grp2idx(categorical(subclass_labels_subset));
        
        % Compute silhouette scores and store median value
        silhouette_scores = silhouette(pc_subset, numeric_labels);
        cross_species_silhouette_matrix(i, j) = median(silhouette_scores);
    end
end

% % Set diagonal to zero to appear white (self-comparisons not meaningful)
% cross_species_silhouette_matrix(eye(size(cross_species_silhouette_matrix)) == 1) = 0;

% Set color scale limits from 0.5 to 1 for consistent interpretation
clims = [0, 1];

% Define white-red-blue colormap
custom_colormap = [linspace(1,1,256)', linspace(1,0,256)', linspace(1,0,256)';  % White to red
                   linspace(1,0,256)', zeros(256,1), linspace(0,1,256)'];       % Red to blue

% Plot heatmap for cross-species subclass comparison
figure;
h_cross_species = heatmap(opossum_subclasses, mouse_subclasses, round(cross_species_silhouette_matrix,2)', ...
    'Colormap', custom_colormap, 'ColorLimits', clims);
title('Cross-Species Mouse vs Opossum IT Subclass Silhouette Heatmap');
xlabel('Opossum Subclass');
ylabel('Mouse Subclass');
h_cross_species.CellLabelFormat = '%.2f';  % Round values to 2 decimal places
h_cross_species.GridVisible = 'off';  % Remove gridlines for a cleaner look

%% Heatmaps (Cross-Species Opossum / Macaque Silhouette Comparison)

% Load data from CSV for both species projected into the same PC space
mouse_data = readtable('E:/macaque_pc_data_shared_with_subclass.csv');
opossum_data = readtable('E:/opossum_pc_data_shared_with_subclass.csv');

% Extract the first three PCs for both species
mouse_pc_matrix = [mouse_data.PC_1, mouse_data.PC_2, mouse_data.PC_3];
opossum_pc_matrix = [opossum_data.pca_1, opossum_data.pca_2, opossum_data.pca_3];

% Convert subclass to categorical and assign numerical labels
mouse_subclasses = unique(mouse_data.subclass);
opossum_subclasses = unique(opossum_data.subclass);

% Combine data for cross-species comparison
combined_pc_matrix = [mouse_pc_matrix; opossum_pc_matrix];
combined_subclass_labels = [mouse_data.subclass; opossum_data.subclass];

% Get unique subclass names across both species
combined_subclasses = unique(combined_subclass_labels);

% Initialize matrix to store silhouette scores for cross-species subclass pairs
cross_species_silhouette_matrix = zeros(length(mouse_subclasses), length(opossum_subclasses));

% Compute silhouette scores for each subclass pair across species
for i = 1:length(mouse_subclasses)
    for j = 1:length(opossum_subclasses)
        % Select cells belonging to the mouse and opossum subclasses
        idx = strcmp(combined_subclass_labels, mouse_subclasses{i}) | strcmp(combined_subclass_labels, opossum_subclasses{j});
        subclass_labels_subset = combined_subclass_labels(idx);
        pc_subset = combined_pc_matrix(idx, :);
        
        % Convert subclass labels to numerical indices
        numeric_labels = grp2idx(categorical(subclass_labels_subset));
        
        % Compute silhouette scores and store median value
        silhouette_scores = silhouette(pc_subset, numeric_labels);
        cross_species_silhouette_matrix(i, j) = median(silhouette_scores);
    end
end

% % Set diagonal to zero to appear white (self-comparisons not meaningful)
% cross_species_silhouette_matrix(eye(size(cross_species_silhouette_matrix)) == 1) = 0;

% Set color scale limits from 0.5 to 1 for consistent interpretation
clims = [0, 1];

% Define white-red-blue colormap
custom_colormap = [linspace(1,1,256)', linspace(1,0,256)', linspace(1,0,256)';  % White to red
                   linspace(1,0,256)', zeros(256,1), linspace(0,1,256)'];       % Red to blue

% Plot heatmap for cross-species subclass comparison
figure;
h_cross_species = heatmap(mouse_subclasses, opossum_subclasses, round(cross_species_silhouette_matrix,2)', ...
    'Colormap', custom_colormap, 'ColorLimits', clims);
title('Cross-Species Mouse vs Opossum IT Subclass Silhouette Heatmap');
xlabel('Opossum Subclass');
ylabel('Mouse Subclass');
h_cross_species.CellLabelFormat = '%.2f';  % Round values to 2 decimal places
h_cross_species.GridVisible = 'off';  % Remove gridlines for a cleaner look

%% Heatmaps (Cross-Species Silhouette Comparison with Permutation Test)

% Load data from CSV for both species projected into the same PC space
mouse_data = readtable('E:/mouse_pc_data_shared_subsample_with_subclass.csv');
opossum_data = readtable('E:/opossum_pc_data_shared_subsample_with_subclass.csv');

% Extract the first three PCs for both species
mouse_pc_matrix = [mouse_data.pca_1, mouse_data.pca_2, mouse_data.pca_3];
opossum_pc_matrix = [opossum_data.pca_1, opossum_data.pca_2, opossum_data.pca_3];

% Convert subclass to categorical and assign numerical labels
mouse_subclasses = unique(mouse_data.subclass);
opossum_subclasses = unique(opossum_data.subclass);

% Combine data for cross-species comparison
combined_pc_matrix = [mouse_pc_matrix; opossum_pc_matrix];
combined_subclass_labels = [mouse_data.subclass; opossum_data.subclass];

% Get unique subclass names across both species
combined_subclasses = unique(combined_subclass_labels);

% Initialize matrices to store silhouette scores and p-values
cross_species_silhouette_matrix = zeros(length(mouse_subclasses), length(opossum_subclasses));
p_value_matrix = zeros(length(mouse_subclasses), length(opossum_subclasses));

% Number of permutations for randomization test
num_permutations = 100;
null_distributions = cell(length(mouse_subclasses), length(opossum_subclasses));

% Compute silhouette scores for each subclass pair across species with permutation testing
for i = 1:length(mouse_subclasses)
    for j = 1:length(opossum_subclasses)
        % Select cells belonging to the mouse and opossum subclasses
        idx = strcmp(combined_subclass_labels, mouse_subclasses{i}) | strcmp(combined_subclass_labels, opossum_subclasses{j});
        subclass_labels_subset = combined_subclass_labels(idx);
        pc_subset = combined_pc_matrix(idx, :);
        
        % Convert subclass labels to numerical indices
        numeric_labels = grp2idx(categorical(subclass_labels_subset));
        
        % Compute observed silhouette score
        observed_silhouette_scores = silhouette(pc_subset, numeric_labels);
        observed_median_silhouette = median(observed_silhouette_scores);
        cross_species_silhouette_matrix(i, j) = observed_median_silhouette;

        % Permutation test: shuffle labels 1000 times
        shuffled_silhouettes = zeros(num_permutations, 1);
        for k = 1:num_permutations
            shuffled_labels = numeric_labels(randperm(length(numeric_labels)));
            shuffled_silhouettes(k) = median(silhouette(pc_subset, shuffled_labels));
        end
        
        % Store null distribution for plotting
        null_distributions{i, j} = shuffled_silhouettes;
        
        % Compute p-value (proportion of shuffled scores greater than or equal to observed)
        p_value_matrix(i, j) = mean(shuffled_silhouettes >= observed_median_silhouette);
    end
end

% Set color scale limits from 0.5 to 1 for consistent interpretation
clims = [0.5, 1];

% Define white-red-blue colormap
custom_colormap = [linspace(1,1,256)', linspace(1,0,256)', linspace(1,0,256)';  % White to red
                   linspace(1,0,256)', zeros(256,1), linspace(0,1,256)'];       % Red to blue

% Plot heatmap for cross-species subclass comparison
figure;
h_cross_species = heatmap(mouse_subclasses, opossum_subclasses, round(cross_species_silhouette_matrix,2), ...
    'Colormap', custom_colormap, 'ColorLimits', clims);
title('Cross-Species Mouse vs Opossum IT Subclass Silhouette Heatmap');
xlabel('Opossum Subclass');
ylabel('Mouse Subclass');
h_cross_species.CellLabelFormat = '%.2f';  % Round values to 2 decimal places
h_cross_species.GridVisible = 'off';  % Remove gridlines for a cleaner look

%% Plot histograms of null distributions with observed values
figure;
tiledlayout(length(mouse_subclasses), length(opossum_subclasses), 'TileSpacing', 'compact');

for i = 1:length(mouse_subclasses)
    for j = 1:length(opossum_subclasses)
        nexttile;
        
        % Get null distribution and observed value
        histogram(null_distributions{i, j}, 'Normalization', 'probability', 'FaceColor', [0.7 0.7 0.7]);
        hold on;
        xline(cross_species_silhouette_matrix(i, j), 'r', 'LineWidth', 2);
        hold off;
        
        % Annotate plot
        title(sprintf('%s vs %s\np=%.3f', mouse_subclasses{i}, opossum_subclasses{j}, p_value_matrix(i, j)));
        xlabel('Silhouette Score');
        ylabel('Frequency');
    end
end
