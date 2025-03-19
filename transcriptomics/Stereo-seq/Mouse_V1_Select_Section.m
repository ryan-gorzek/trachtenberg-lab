
% Load and plot original data
y_coords = table2array(readtable("E:/STOmics/seurat/Mouse/raw/x_coords.csv")) * -1;
x_coords = table2array(readtable("E:/STOmics/seurat/Mouse/raw/y_coords.csv")) * -1;
barcodes = string(table2array(readtable("E:/STOmics/seurat/Mouse/raw/barcodes.tsv", "FileType", "text", 'Delimiter', '\t')));

% Plot the spatial transcriptomics data
figure;
scatter(x_coords, y_coords, 1, 'filled'); % Adjust marker size if needed
xlabel('X Coordinate');
ylabel('Y Coordinate');
title('Select an ROI using ginput');
hold on;

% Interactive ROI selection
disp('Click to define an ROI. Press Enter when done.');

x_roi = [];
y_roi = [];

while true
    [x, y, button] = ginput(1); % Get one point at a time
    
    if isempty(button) % If Enter is pressed, exit loop
        break;
    end
    
    % Store the clicked point
    x_roi = [x_roi; x];
    y_roi = [y_roi; y];
    
    % Plot the clicked points
    plot(x_roi, y_roi, 'ro-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
end

% Close the polygon by connecting last to first point
if length(x_roi) > 2
    x_roi = [x_roi; x_roi(1)];
    y_roi = [y_roi; y_roi(1)];
    plot(x_roi, y_roi, 'r-', 'LineWidth', 2);
end

% Find points inside the ROI
in_roi = inpolygon(x_coords, y_coords, x_roi, y_roi);

% Extract barcodes of selected cells
selected_barcodes = barcodes(in_roi);

% Save selected barcodes to a file
writematrix(selected_barcodes, "E:/STOmics/seurat/Mouse/selected_barcodes.csv");

disp('Selected barcodes saved to selected_barcodes.txt');
