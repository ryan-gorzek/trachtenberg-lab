
% Load and plot original data
y_coords = table2array(readtable("E:/STOmics/seurat/Opossum/raw/x_coords.csv")) * -1;
x_coords = table2array(readtable("E:/STOmics/seurat/Opossum/raw/y_coords.csv")) * -1;

% Initialize parameters array
region_params = [];

% Example region parameters (manually defined or in a loop)
angle_deg = 7; % rotation angle in degrees
theta = deg2rad(angle_deg);

% Rotate coordinates
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
rotated_coords = R * [x_coords'; y_coords'];
x_rot = rotated_coords(1, :)';
y_rot = rotated_coords(2, :)';

% Plot rotated data for region estimation
figure; set(gcf, "Color", "w");
scatter(x_rot, y_rot, 5, '.r');
title(sprintf('Rotated by %d°', angle_deg));
set(gca, "Box", "on", "LineWidth", 2, "TickDir", "out", "FontSize", 10);

% 🔻 Manually define region boundaries after visual inspection:
x1 = -7000; x2 = -5000;
y1 = -14600; y2 = -13050;

% Save parameters
region_params = [region_params; angle_deg, x1, x2, y1, y2];

%% Plot original data with region overlay
figure; set(gcf, "Color", "w");
scatter(x_coords, y_coords, 5, '.r');
title('Original Coordinates with Regions');
set(gca, "Box", "on", "LineWidth", 2, "TickDir", "out", "FontSize", 10);
hold on;

% Loop over regions to plot them back on the original (unrotated) plot
for i = 1:size(region_params, 1)
    angle_deg = region_params(i, 1);
    theta = deg2rad(-angle_deg); % reverse rotation

    % Define rectangle in rotated space
    rect_x = [region_params(i,2), region_params(i,3), region_params(i,3), region_params(i,2), region_params(i,2)];
    rect_y = [region_params(i,4), region_params(i,4), region_params(i,5), region_params(i,5), region_params(i,4)];

    % Rotate rectangle back to original space
    R_inv = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    rect_coords = R_inv * [rect_x; rect_y];

    % Plot rectangle
    plot(rect_coords(1,:), rect_coords(2,:), '-k', 'LineWidth', 1.5);
end

hold off;
