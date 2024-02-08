
quad = load_quad_data("pv_vip_00", "401", "000");

%%

figure; tiledlayout(4, 6);

raw = movmean(double(quad), 15);
drv = abs([0 diff(raw)]);

nexttile([1, 3]); plot(raw); axis tight;
nexttile([1, 3]); plot(drv); axis tight;
nexttile([3, 3]); histogram(corr(horzcat(dFF.cells.data), raw'));
nexttile([3, 3]); histogram(corr(horzcat(dFF.cells.data), drv'));
