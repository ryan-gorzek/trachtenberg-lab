
%%

[eye, imgs] = load_eye_data("pv_vip_00", "310", "000");

%% make a movie

v = VideoWriter("eye_movie_as", "Archival");
open(v);
for fr = 1:size(imgs, 3)
    writeVideo(v, imgs(:, :, fr));
end
close(v);

%%

% area = horzcat(eye.Area);
% % area = area - median(area);
% % area = area(150 : end - 150);
% area = medfilt1(area, 200, "truncate");
% quad = abs(diff(double(quad_data)));
% quad = quad./max(quad);
% figure; hold on; plot(quad*2); plot(zscore(area(1:end-1))); 

%%

centroid = vertcat(eye.Centroid);

figure;
for fr = 1:size(centroid, 1)
imshow(imgs(:, :, fr)); hold on;
scatter(centroid(fr, 1), centroid(fr, 2), "+r");
drawnow;
end
% th = 0:pi/50:2*pi;
% xunit = r * cos(th) + x;
% yunit = r * sin(th) + y;
% h = plot(xunit, yunit);

%%

figure;
imshow(data_crop); hold on;
[x, y] = ginput(1);
scatter(x, y, "+r");
