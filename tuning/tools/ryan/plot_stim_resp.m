function plot_stim_resp(pop, frameon, window)

nstim = numel(frameon);
ncells = size(pop, 2);
nframes = numel(window(1) : window(2));
resp = zeros(ncells, nframes, nstim);
for n = 1:nstim
    cw = frameon(n) + window(1) : frameon(n) + window(2);
    resp(:, :, n) = pop(cw, :)';
end
if ncells == 1
    resp = permute(resp, [3, 2, 1]);
end
xline(-15, "--k", "LineWidth", 2);
xline(0, "--k", "LineWidth", 2);
xline(15, "--k", "LineWidth", 2);
plot(window(1) : window(2), mean(resp, 3), "LineWidth", 1, "Color", [0.8, 0.8, 0.8, 0.7]);
% plot(window(1) : window(2), mean(resp, [1, 3]), "LineWidth", 4, "Color", "k");
plot(window(1) : window(2), std(squeeze(resp)), "LineWidth", 4, "Color", "k");
axis tight;
xlabel("Frames from Stimulus Onset");
set(gca, "FontSize", 9, ...
         "XColor", "k", "YColor", "k", ...
         "TickDir", "out", "TickLength", [0.01, 0.01], ...
         "box", "off", ...
         "LineWidth", 1);
end