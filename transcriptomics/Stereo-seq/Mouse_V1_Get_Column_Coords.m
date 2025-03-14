
y_coords = table2array(readtable("E:/STOmics/seurat/Mouse/raw/x_coords.csv")) * -1;
x_coords = table2array(readtable("E:/STOmics/seurat/Mouse/raw/y_coords.csv")) * -1;

figure; set(gcf, "Color", "w");
s = scatter(x_coords, y_coords, 5, '.r');
set(gca, "Box", "on", "LineWidth", 2, "TickDir", "out", "FontSize", 10);
pb = pbaspect;
set(gcf, "Position", [813, 146.1, 1360.6, 888.6]);
pbaspect(pb);

