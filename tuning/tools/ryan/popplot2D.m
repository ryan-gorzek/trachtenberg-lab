function f = popplot2D(dFF, spikes)

f = figure; tiledlayout(2, 2, "TileSpacing", "tight");

peakvals_dFF = vertcat(dFF.kernels.peakval);
for d = 1:dFF.ndims
    if dFF.dimscales(d) == "log"
        peakvals_dFF(:, d) = log10(peakvals_dFF(:, d));
    end
end

nexttile; histogram(peakvals_dFF(:, 1), spikes.viz.binedges{1});
xticks(spikes.viz.ticks.cont{1}); xticklabels(spikes.viz.ticklabels{1});
title(dFF.dimnames{1});
ylabel(["dFF"; "# of Cells"]);
xline(median(peakvals_dFF(:, 1)), "--r", "LineWidth", 2);
setStyle();
nexttile; histogram(peakvals_dFF(:, 2), dFF.viz.binedges{2});
xticks(dFF.viz.ticks.cont{2}); xticklabels(dFF.viz.ticklabels{2});
title(dFF.dimnames{2});
xline(median(peakvals_dFF(:, 2)), "--r", "LineWidth", 2);
setStyle();

peakvals_spikes = vertcat(spikes.kernels.peakval);
for d = 1:spikes.ndims
    if spikes.dimscales(d) == "log"
        peakvals_spikes(:, d) = log10(peakvals_spikes(:, d));
    end
end

nexttile; histogram(peakvals_spikes(:, 1), spikes.viz.binedges{1});
xticks(spikes.viz.ticks.cont{1}); xticklabels(spikes.viz.ticklabels{1});
ylabel(["Spikes"; "# of Cells"]);
xline(median(peakvals_dFF(:, 1)), "--r", "LineWidth", 2);
setStyle();
nexttile; histogram(peakvals_spikes(:, 2), spikes.viz.binedges{2});
xticks(spikes.viz.ticks.cont{2}); xticklabels(spikes.viz.ticklabels{2});
xline(median(peakvals_spikes(:, 2)), "--r", "LineWidth", 2);
setStyle();

end
