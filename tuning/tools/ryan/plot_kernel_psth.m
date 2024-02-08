function plot_kernel_psth(exp_psths, dim)

kernel_psth = exp_psths.kern;
[nrow, ncol] = size(kernel_psth);
figure; tiledlayout(nrow, ncol, "TileSpacing", "tight");
kernel_psth = reshape(kernel_psth', 1, []);
for p = 1:numel(kernel_psth)
    resp = kernel_psth{p};
    nexttile; hold on;
    xline(-15, "--k", "LineWidth", 2);
    xline(0, "--k", "LineWidth", 2);
    xline(15, "--k", "LineWidth", 2);
    plot(exp_psths.window(1) : exp_psths.window(2), resp, "LineWidth", 1, "Color", [0.8, 0.8, 0.8, 0.7]);
    plot(exp_psths.window(1) : exp_psths.window(2), mean(resp, 1), "LineWidth", 4, "Color", "k");
    axis off;
end

end