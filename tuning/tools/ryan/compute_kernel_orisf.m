function [kernel_out_raw, kernel_out_filt, stats] = compute_kernel_orisf(datamatrix, frameon, frameoff, params)
% datamatrix is 

ncells = size(datamatrix, 2);
nstims = numel(frameon);
ntaus = 20;

r = zeros([size(params.pp1), ntaus, ncells]);
N = zeros(size(params.pp1));

for s = 1:nstims
    r(params.stimidx(s, 2), params.stimidx(s, 1), params.stimidx(s, 3), :, :) = ...
        squeeze(r(params.stimidx(s, 2), params.stimidx(s, 1), params.stimidx(s, 3), :, :)) + ...
        datamatrix(frameon(s) - 2 : frameon(s) + ntaus - 3, :);
    N(params.stimidx(s, 2), params.stimidx(s, 1), params.stimidx(s, 3)) = ...
        N(params.stimidx(s, 2), params.stimidx(s, 1), params.stimidx(s, 3)) + 1;
end

for t = 1:ntaus
    for n = 1:ncells
        r(:, :, :, t, n) = r(:, :, :, t, n) ./ N;
    end
end

%%%% complexity (F1/F0), cell-wise in loop below
R = r;
R(isnan(R)) = 0;
F1F0 = abs(fft(R, [], 3));
F1F0 = squeeze(F1F0(:, :, 2, :, :) ./ F1F0(:, :, 1, :, :));

r = squeeze(mean(r, 3, "omitnan")); % average across spatial phases
rfilt = nan(size(r));

h = fspecial('gauss', 5, 1);
k = 0;
for t = 1:ntaus
    for n = 1:ncells
        % extend kernel by replicating edges for filtering
        rf = squeeze(r(:, :, t, n));
        rf2 = [rf(end-1, :); rf(end, :); rf; rf(1, :); rf(2, :)];
        rf3 = [rf2(:, 1), rf2(:, 1), rf2, rf2(:, end), rf2(:, end)];
        rf4 = filter2(h, rf3, 'valid');
        rfilt(:, :, t, n) = rf4;
        k = k + 1;
    end
end

kernel_out_raw = nan([size(params.pp1, [1, 2]), ncells]);
kernel_out_filt = nan([size(params.pp1, [1, 2]), ncells]);
stats = struct("SNR", cell(ncells, 1), ...
               "tmax", cell(ncells, 1), ...
               "F1F0", cell(ncells, 1));

SNR = nan(ncells, 1);
for i = 1:ncells
    %%%% compute kernels
    zraw = squeeze(r(:, :, :, i));
    zfilt = squeeze(rfilt(:, :, :, i));
    qfilt = reshape(zfilt, params.np1 * params.np2, []);
    kfilt = std(qfilt);
    [~, tr] = max(kfilt(7:11));
    kernel_out_raw(:, :, i) = squeeze(zraw(:, :, tr + 6));
    kernel_out_filt(:, :, i) = squeeze(zfilt(:, :, tr + 6));
    %%%% compute stats
    [kmax, tmax] = max(kfilt);
    if kmax < 0.005
        snr = mean(kfilt(8:10)) / mean(kfilt);
        response = mean(kfilt(8:10));
        noise = mean(kfilt);
        noise_std = std(kfilt);
    else
        snr = mean(kfilt(8:10)) / mean(kfilt([1 2 3 16:20]));
        SNR(i) = snr;
        response = mean(kfilt(8:10));
        noise = mean(kfilt([1 2 3 16:20]));
        noise_std = std(kfilt([1 2 3 16:20]));
    end
    stats(i).SNR = snr;
    stats(i).tmax = tmax;
    % F1/F0
    [idx_ori, idx_sf] = find(kernel_out_filt(:, :, i) == max(kernel_out_filt(:, :, i), [], "all"));
    stats(i).F1F0 = 2 * F1F0(idx_ori, idx_sf, tmax, i);
    if any(size(stats(i).F1F0) > 1), stats(i).F1F0 = NaN; end

%     stat(i).q = q;
%     stat(i).k = k;
%     stat(i).response = response;
%     stat(i).noise = noise;
%     stat(i).noise_std = noise_std;
%     stat(i).kmax = kmax;

end

% assignin("base", "stats", stats); assignin("base", "rfilt", rfilt); pause;

end