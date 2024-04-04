
function eye = sbxeyemotion(fn, rad_range, concat_fpath)

fprintf("Processing pupil data...\n");

fn = char(fn);
load(fn, 'data');                  % should be a '*_eye.mat' file

if size(rad_range, 1) == 1
    % Edit patrick - allow rad_range to change as a function of time for
    % more robust estimation, esp. with whiskers
    rad_range = ones(size(data, 4), 1) * rad_range;
end

data = squeeze(data);      % the raw images...

% check for center coordinates saved to suite2p plane0 folder
f_eye_center = strcat(concat_fpath, filesep, "suite2p", filesep, "plane0", filesep, "eye_center.csv");
if isfile(f_eye_center)
    eye_center = readmatrix(f_eye_center);
    xc = eye_center(1); yc = eye_center(2);
else
    f = figure; imshow(data(:, :, 1)); hold on;
    [xc, yc] = ginput(1);
    scatter(xc, yc, "+r");
    pause(2);
    close(f);
end

W = 50;

x_start = floor(xc - W); x_end = ceil(xc + W);
y_start = floor(yc - W); y_end = ceil(yc + W);

if any([x_start, y_start] < 1) || (x_end > size(data, 2)) || (y_end > size(data, 1))
    error("Cannot crop image, choose a new center...");
else
    writematrix([xc, yc], f_eye_center);
    data_crop = data(y_start : y_end, x_start : x_end, :);
end

% Edit Patrick: changed this so that struct doesn't change size every time
% while maintaining backwards compatibility

A = cell(size(data_crop, 3), 1);
B = cell(size(data_crop, 3), 1);
C = cell(size(data_crop, 3), 1);
D = cell(size(data_crop, 3), 1);
for n = 1:size(data_crop, 3)
    A{n} = [0, 0];
    B{n} = 0;
    C{n} = 0;
    D{n} = {};
end

eye = struct('Centroid', A, 'Radius', B, 'Area', C, 'Data', D);
parfor n = 1:size(data, 3)
    warning off;
    [center, radii, metric] = imfindcircles(squeeze(data_crop(:, :, n)), rad_range(n, :), 'Sensitivity', 0.95);
    if isempty(center)
        eye(n).Centroid = [NaN NaN];    % could not find anything...
        eye(n).Radius = NaN;
        eye(n).Area = NaN;
    else
        [~, idx] = sort(metric, "descend");          % pick the circle with best score
        idx = idx(1:min(2, numel(idx))); % choose top 2
        dist_xy = center(idx, :) - size(data_crop, [1, 2]) ./ 2;
        dist_euc = sqrt(sum(dist_xy .^ 2, 2));
        if nnz(dist_euc < 25) == 1 % if there's only one circle close to the center, get that one
            [~, idx_get] = min(dist_euc);
        else % if there's two, choose the higher score
            idx_get = idx(1);
        end
        eye(n).Centroid = center(idx(idx_get), :) + [xc - W, yc - W]; % account for cropping
        eye(n).Radius = radii(idx(idx_get));
        eye(n).Area = pi*radii(idx(idx_get))^2;
        eye(n).Data = {center, radii, metric};
    end
end

save(fn, 'eye', '-append');     % append the motion estimate data...
