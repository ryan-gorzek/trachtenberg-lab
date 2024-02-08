
[eye, imgs] = load_eye_data("pv_vip_00", "401", "000");

%%

figure; hold on;
for fr = 18775
% for fr = 4143
imshow(imgs(:, :, fr));

center = eye(fr).Centroid;
radius = sqrt(eye(fr).Area / pi);
viscircles(center, radius);

drawnow;

end

%%

figure; tiledlayout(4, 6);

filt = @(x) medfilt1(x(1:end-2), 15, "truncate");
% filt = @(x) movmean(x(1:end-2), 15);
diam = filt(2 .* vertcat(eye.Radius));
area = filt(vertcat(eye.Area));

nexttile([1, 3]); plot(diam); axis tight;
nexttile([1, 3]); plot(area); axis tight;
nexttile([3, 3]); histogram(corr(horzcat(dFF.cells.data), diam));
nexttile([3, 3]); histogram(corr(horzcat(dFF.cells.data), area));

%%

figure; tiledlayout(4, 6);

raw = medfilt1(dFF.eye.data(2:end-1), 15, "truncate");
drv = abs([0; diff(raw)]);

nexttile([1, 3]); plot(raw); axis tight;
nexttile([1, 3]); plot(drv); axis tight;
nexttile([3, 3]); histogram(corr(horzcat(dFF.cells.data), raw));
nexttile([3, 3]); histogram(corr(horzcat(dFF.cells.data), drv));

%%

figure; t = tiledlayout(2, 1);

cellnum = 704;
dff = circshift(dFF.cells(cellnum).data, 5);

nexttile; plot(dff); axis tight;
nexttile; plot(raw); axis tight;

title(t, sprintf("Corr = %.3f", corr(dff, raw)));
