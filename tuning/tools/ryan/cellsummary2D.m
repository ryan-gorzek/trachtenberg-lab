function f = cellplot2D(expstruct)


f = figure; set(gcf, "Color", "w"); tiledlayout(2, 6, "TileSpacing", "tight");

%%%% dFF
% trace
nexttile([1, 2]); plot(expstruct.cells(id).data, "LineWidth", 0.5, "Color", "k"); title("\DeltaF / F"); axis tight;
% PSTH
nexttile; hold on; plot_stim_resp(expstruct.cells(id).data, frame_on, dFF.viz.window);
% kernel
nexttile; imagesc(dFF.kernels(id).kern); colorbar;
xticks(expstruct.viz.ticks.disc{2}); xticklabels(expstruct.viz.ticklabels{2}); xlabel(expstruct.viz.labels{2});
yticks(expstruct.viz.ticks.disc{1}); yticklabels(expstruct.viz.ticklabels{1}); ylabel(expstruct.viz.labels{1});
% dim1 curve
nexttile; plot(expstruct.viz.ticks.cont{1}, expstruct.kernels(id).curve{1}, "-ok"); axis tight;
xticks(expstruct.viz.ticks.cont{1}); xticklabels(expstruct.viz.ticklabels{1}); xlabel(expstruct.viz.labels{1});
% dim2 curve
nexttile; plot(expstruct.viz.ticks.cont{2}, expstruct.kernels(id).curve{2}, "-ok"); axis tight;
xticks(expstruct.viz.ticks.cont{2}); xticklabels(expstruct.viz.ticklabels{2}); xlabel(expstruct.viz.labels{2});

end