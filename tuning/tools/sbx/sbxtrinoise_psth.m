function r = sbxtrinoise_psth(fname, info)

% analyze tri noise experiment

% edit Luis 1/10/17. Allow rigid and nonrigid .signals files as input..
    % -----
    if contains(fname,'rigid') % search for rigid in filename
        si = strfind(fname,'_'); 
        fnamelog = fname( 1:si(end)-1); % remove it
    else 
        fnamelog = fname;
    end
    % -----

%%
log = sbxreadtrinoiselog(fnamelog); % read log

load([fname, '_suite.signals'],'-mat', 'spks');    % load signals
if ~ismatrix(spks)
    spks = squeeze(spks(1,:,:));     % keep green channel
    % sig =  squeeze(sig(1,:,:));
end

dsig = spks;

ncell = size(dsig, 2);
nstim = length(log);

[nrow, ncol] = size(log{1}{2});

X = zeros(nstim,nrow,ncol); % data matrix
for i = 1:nstim
    X(i,:,:) = log{i}{2};
end

frame = info.evt.photodiode_off.frame;

ntau = 20;

kern_on = zeros(ncell,ntau,nrow,ncol);
kern_off = zeros(ncell,ntau,nrow,ncol);

stim_on  = zeros(nrow,ncol);
stim_off = zeros(nrow,ncol);

for i = 1:nstim
    stim_on = stim_on   + double(squeeze(X(i,:,:)) >  0);
    stim_off = stim_off + double(squeeze(X(i,:,:)) <  0);
end

parfor n=1:ncell
    for i=1:nstim
        for t=1:ntau
            kern_on(n,t,:,:)  = squeeze(kern_on(n,t,:,:))   + dsig(frame(i)+t-3,n)   * double(squeeze(X(i,:,:)) >  0);
            kern_off(n,t,:,:) = squeeze(kern_off(n,t,:,:))  + dsig(frame(i)+t-3,n)   * double(squeeze(X(i,:,:)) <  0);
        end
    end
end

% normalize

for n=1:ncell
    for t = 1:ntau
        kern_on(n,t,:,:) = squeeze(kern_on(n,t,:,:))./stim_on;
        kern_off(n,t,:,:) = squeeze(kern_off(n,t,:,:))./stim_off;
    end
end

% kernel analysis

r = struct([]);
H = fspecial('gauss',3,.6);
for n = 1:ncell

    fprintf("Fitting trinoise kernel for cell #%i out of %i...\n", n, ncell);
    tic;

    k = squeeze(kern_on(n,:,:,:));
    k = reshape(k,ntau,[]);
    s = sqrt(sum(k.^2,2));
    snr_on = s/s(1);
    [snr_on_max,tmax_on] = max(snr_on);
    k_on  = filter2(H,squeeze(kern_on(n,tmax_on,:,:)),'valid');
    [xx,yy] = meshgrid(1:size(k_on,2),1:size(k_on,1));
    [fitres,zfit] = sbxgaussfit(xx,yy,k_on,1);
    
    r(n).snr_on = snr_on;
    r(n).snr_on_max = snr_on_max;
    r(n).tmax_on = tmax_on;
    r(n).kern_on = k_on;
    r(n).kurt_on = kurtosis(k_on(:));
    r(n).kern_on_fit = zfit;
    r(n).kern_on_param = fitres;
    rho = corrcoef(zfit(:),k_on(:));
    r(n).cc_on = rho(1,2);

    [xx,yy] = meshgrid(1:size(k_on,2),1:size(k_on,1));
    k_on = (k_on - max(k_on(:))/sqrt(2));
    k_on = k_on .* (k_on > 0);
    r(n).on_x = sum(k_on(:).*xx(:))/sum(k_on(:));
    r(n).on_y = sum(k_on(:).*yy(:))/sum(k_on(:));
    r(n).on_amp = max(k_on(:));

    [xx,yy] = meshgrid(1:size(r(n).kern_on_fit,2),1:size(r(n).kern_on_fit,1));
    k_on_fit = (r(n).kern_on_fit - max(r(n).kern_on_fit(:))/sqrt(2));
    k_on_fit = k_on_fit .* (k_on_fit > 0);
    r(n).on_x_fit = sum(k_on_fit(:).*xx(:))/sum(k_on_fit(:));
    r(n).on_y_fit = sum(k_on_fit(:).*yy(:))/sum(k_on_fit(:));
    r(n).on_amp_fit = max(k_on_fit(:));
    
    k = squeeze(kern_off(n,:,:,:));
    k = reshape(k,ntau,[]);
    s = sqrt(sum(k.^2,2));
    snr_off = s/s(1);
    [snr_off_max,tmax_off] = max(snr_off);
    k_off  = filter2(H,squeeze(kern_off(n,tmax_off,:,:)),'valid');
    [xx,yy] = meshgrid(1:size(k_on,2),1:size(k_on,1));
    [fitres,zfit] = sbxgaussfit(xx,yy,k_off,1);
    
    r(n).snr_off = snr_off;
    r(n).snr_off_max = snr_off_max;
    r(n).tmax_off = tmax_off;
    r(n).kern_off = k_off;
    r(n).kurt_off = kurtosis(k_off(:));
    r(n).kern_off_fit = zfit;
    r(n).kern_off_param = fitres;
    rho = corrcoef(zfit(:),k_off(:));
    r(n).cc_off = rho(1,2);

    [xx,yy] = meshgrid(1:size(k_off,2),1:size(k_off,1));
    k_off = (k_off - max(k_off(:))/sqrt(2));
    k_off = k_off .* (k_off > 0);
    r(n).off_x = sum(k_off(:).*xx(:))/sum(k_off(:));
    r(n).off_y = sum(k_off(:).*yy(:))/sum(k_off(:));
    r(n).off_amp = max(k_off(:));

    [xx,yy] = meshgrid(1:size(r(n).kern_off_fit,2),1:size(r(n).kern_off_fit,1));
    k_off_fit = (r(n).kern_off_fit - max(r(n).kern_off_fit(:))/sqrt(2));
    k_off_fit = k_off_fit .* (k_off_fit > 0);
    r(n).off_x_fit = sum(k_off_fit(:).*xx(:))/sum(k_off_fit(:));
    r(n).off_y_fit = sum(k_off_fit(:).*yy(:))/sum(k_off_fit(:));
    r(n).off_amp_fit = max(k_off_fit(:));

    toc;

end

save([fname, '.trinoise_kern'], 'r');
