function f = cellplot2D(dFF, spikes, id)

if isfield(dFF, "psths"), ncol = 6; else, ncol = 5; end
f = figure; set(gcf, "Color", "w"); tiledlayout(2, ncol, "TileSpacing", "tight");

%%%% dFF
% trace
nexttile([1, 2]); plot(dFF.cells(id).data, "LineWidth", 0.5, "Color", "k"); title("\DeltaF / F"); axis tight;
% PSTH
if isfield(dFF, "psths")
    nexttile; hold on; plot_stim_resp(dFF.cells(id).data, vertcat(dFF.stims.onset), dFF.viz.window);
end
% kernel
nexttile; imagesc(dFF.kernels(id).kernsmooth); pbaspect([flip(size(dFF.kernels(id).kern)), 1]); colorbar;
xticks(dFF.viz.ticks.disc{2}); xticklabels(dFF.viz.ticklabels{2}); xlabel(dFF.viz.labels{2});
yticks(dFF.viz.ticks.disc{1}); yticklabels(dFF.viz.ticklabels{1}); ylabel(dFF.viz.labels{1});
% dim1 curve
nexttile; plot(dFF.viz.ticks.cont{1}, dFF.kernels(id).curve{1}, "-ok"); axis tight;
xticks(dFF.viz.ticks.cont{1}); xticklabels(dFF.viz.ticklabels{1}); xlabel(dFF.viz.labels{1});
% dim2 curve
nexttile; plot(dFF.viz.ticks.cont{2}, dFF.kernels(id).curve{2}, "-ok"); axis tight;
xticks(dFF.viz.ticks.cont{2}); xticklabels(dFF.viz.ticklabels{2}); xlabel(dFF.viz.labels{2});

%%%% spikes
% trace
nexttile([1, 2]); plot(spikes.cells(id).data, "LineWidth", 0.5, "Color", "k"); title("Inferred Spikes"); axis tight;
% PSTH
if isfield(spikes, "psths")
    nexttile; hold on; plot_stim_resp(spikes.cells(id).data, vertcat(dFF.stims.onset), spikes.viz.window);
end
% kernel
nexttile; imagesc(spikes.kernels(id).kernsmooth); pbaspect([flip(size(spikes.kernels(id).kern)), 1]); colorbar;
xticks(spikes.viz.ticks.disc{2}); xticklabels(spikes.viz.ticklabels{2}); xlabel(spikes.viz.labels{2});
yticks(spikes.viz.ticks.disc{1}); yticklabels(spikes.viz.ticklabels{1}); ylabel(spikes.viz.labels{1});
% dim1 curve
nexttile; plot(spikes.viz.ticks.cont{1}, spikes.kernels(id).curve{1}, "-ok"); axis tight;
xticks(spikes.viz.ticks.cont{1}); xticklabels(spikes.viz.ticklabels{1}); xlabel(spikes.viz.labels{1});
% dim2 curve
nexttile; plot(spikes.viz.ticks.cont{2}, spikes.kernels(id).curve{2}, "-ok"); axis tight;
xticks(spikes.viz.ticks.cont{2}); xticklabels(spikes.viz.ticklabels{2}); xlabel(spikes.viz.labels{2});

end