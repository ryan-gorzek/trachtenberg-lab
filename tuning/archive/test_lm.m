
%%%% add tools to path
addpath(genpath(fullfile("E:/Imaging/code/tools/")));
% addpath(genpath(fullfile("/Users/ryan.gorzek/Imaging/code/tools/")));

subj = "pv_vip_01"; sess = "340"; run = "000"; stim = "battery4";

% load calcium data
[rtdata, stdata, xy, grid] = load_calcium_data(subj, sess, run);

% load stimulus data
[frame_on, frame_off, params] = load_stimulus_data(subj, sess, run);
if sess == "340", frame_on = frame_on(1:end-1); end
% load quadrature data
quad = load_quad_data(subj, sess, run);

% load eye data
eye = load_eye_data(subj, sess, run);

% create experiment structure
[spikes, dFF] = expstruct(stim, rtdata, stdata, xy, grid, frame_on, frame_off, params, eye, quad);

if stim == "randorisf", SNR = vertcat(spikes.stats.SNR) > spikes.calc.SNR_thr; end

%%

fit = zeros(61, 1);
shifts = -30:30;

for s = 1:size(fit, 1)

    X = circshift(zscore([dFF.regressors, dFF.eye.data, dFF.quad.data]), shifts(s), 1);
    Y = dFF.cells(284).data;

    mdl = fitglm(X, Y);

    fit(s) = mdl.Rsquared.Ordinary;

end

%%

figure; tiledlayout(1, 2);
nexttile; hold on; plot(Y); plot(mdl.Fitted.Response);
nexttile; stem(mdl.Coefficients.Estimate);
