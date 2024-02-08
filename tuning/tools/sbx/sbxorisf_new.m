function stat = sbxorisf_new(fname)

% analyze sparse noise experiment

% edit Luis 1/10/17. Allow rigid and nonrigid .signals files as input..
% edit Liming 1/8/2019.
% edit Liming 6/19/2019 to add f1f0 calculation.
    % -----
    if contains(fname,'rigid') % search for rigid in filename
        si = strfind(fname,'_'); 
        fnamelog = fname( 1:si(end)-1); % remove it
    elseif contains(fname,'suite') % search for suite in filename
        si = strfind(fname,'_'); 
        fnamelog = fname( 1:si(end)-1); % remove it
    else 
        fnamelog = fname;
    end
    % -----

%%
log = sbxreadorisflog(fnamelog); % read log
log = log{1};       % assumes only 1 trial
max_ori = max(abs(log.ori))+1;
max_sphase = max(abs(log.sphase))+1;
max_sper = max(abs(log.sper))+1;

load([fname, '.signals'],'-mat');    % load signals
if(ndims(spks)>2)
    spks = squeeze(spks(1,:,:));
    sig =  squeeze(sig(1,:,:));
end


dsig = spks;

ncell = size(dsig,2);
nstim = size(log,1);

ntau = 20;

r = zeros(max_ori,max_sphase,max_sper,ntau,ncell);
N = zeros(max_ori,max_sphase,max_sper);

disp('Processing...');
for(i=1:nstim)
        r(log.ori(i)+1,log.sphase(i)+1,log.sper(i)+1,:,:) =  ... 
            squeeze(r(log.ori(i)+1,log.sphase(i)+1,log.sper(i)+1,:,:)) + ...
            dsig(log.sbxframe(i)-2:log.sbxframe(i)+ntau-3,:);
         N(log.ori(i)+1,log.sphase(i)+1,log.sper(i)+1) = ...
            N(log.ori(i)+1,log.sphase(i)+1,log.sper(i)+1) + 1; 
end

% Normalize by stim count
for(t=1:ntau)
    for(n = 1:ncell)
        r(:,:,:,t,n) = r(:,:,:,t,n)./N;
    end
end


% average across spatial phase

R = r;
r = squeeze(nanmean(r,2));

% f1/f0

R(isnan(R))=0;
f1f0 = abs(fft(R,[],2));
f1f0 = squeeze(f1f0(:,2,:,:,:) ./ f1f0(:,1,:,:,:));

% now r = r(ori,sper,time,cell)

h = fspecial('gauss',5,1);
disp('Filtering');
k = 0;
for(t=1:ntau)
    for(n = 1:ncell)
        % appears to extend kernel by replicating edges for filtering?
        rf = squeeze(r(:,:,t,n));
        rf2 = [rf(end-1,:); rf(end,:); rf; rf(1,:); rf(2,:)];
        rf3 = [rf2(:,1) rf2(:,1) rf2 rf2(:,end) rf2(:,end)];
%         figure; tiledlayout(1, 4); hold on;
%         nexttile; imagesc(rf);
%         nexttile; imagesc(rf2);
%         nexttile; imagesc(rf3);
        rf4 = filter2(h, rf3, 'valid');
%         nexttile; imagesc(rf4);
%         pause;
        r(:,:,t,n) = rf4;
        k = k+1;
    end
end


% param values

ori = 0:20:170;
sper = 1920./(1.31.^(0:11));   % spatial period in pixels %1.395 gets to .32 cyc/deg
sper_pix = sper;
sper = sper / 15.25;           % 15.25 pixels per degree 
sf = 1./sper;                  % cycles/deg
sphase = linspace(0,360,9);
sphase = sphase(1:end-1);




disp('Statistics...')

% [xx,yy] = meshgrid(-12:12,-12:12);
% zz = xx+1i*yy;
% zz = abs(zz).*exp(1i*angle(zz)*2);

% calculate for each case...

disp(size(r));
for i = 1:size(r,4)
    z = squeeze(r(:,:,:,i));
    q = reshape(z,max_ori*max_sper,[]);
    k = std(q);
%     disp(size(q));
%     figure; plot(mean(q, 1)); hold on; plot(mean(k, 1));
%     pause;
    [kmax,t] = max(k);
    %if ((t<11) && (t>6))
    %snr = kmax/mean(k([1 2 3 16:20]));
    %elseif kmax<0.005
    if kmax<0.005
    snr = mean(k([8:10]))/mean(k);
    response = mean(k([8:10]));
    noise = mean(k);
    noise_std = std(k);
    else    
    snr = mean(k([8:10]))/mean(k([1 2 3 16:20]));
    response = mean(k([8:10]));
    noise = mean(k([1 2 3 16:20]));
    noise_std = std(k([1 2 3 16:20]));
    end
    tmax = t;    
    stat(i).q = q;
    stat(i).k = k;
    stat(i).tmax = tmax;
    stat(i).snr = snr;
    stat(i).response = response;
    stat(i).noise = noise;
    stat(i).noise_std = noise_std;
    stat(i).kmax = kmax;
    %if ((t<12) && (t>6))
    %stat(i).kern = squeeze(z(:,:,t));
    %else
    [kr,tr] = max(k([7:11]));
    stat(i).kern = squeeze(z(:,:,(tr+6)));
    %end

    imagesc(log10(sf),ori,stat(i).kern);
    xlabel('Spatial Frequency (cycles/deg)');
    ylabel('Orientation (deg)');
    xval = get(gca,'xtick');
    l = cell(length(xval),1);
    for k = 1:length(l)
        l{k} = sprintf('%.2f',10^xval(k));
    end
    set(gca,'xticklabel',l);
 
% simple complex dario's style
    
    [idx_ori,idx_sf] = find(stat(i).kern == max(stat(i).kern(:)));
    stat(i).f1f0 = 2*f1f0(idx_ori,idx_sf,tmax,i);
    
% separability measure
    [u,s,v] = svd(stat(i).kern);
    stat(i).lambda12 = s(1,1)/s(2,2);
    
    % energy measure
    
    stat(i).delta = max(stat(i).k)/stat(i).k(1);
    stat(i).sig = stat(i).delta>1.75;
        
    % estimate preferred ori and sf
    
    q = s(1,1);
    s = zeros(size(s));
    s(1,1) = q;
    kest = u*s*v';
    stat(i).kest = kest;
    
%     resp_sf = mean(kest);
%     resp_ori = mean(kest');
    
    [ii,jj] = find(stat(i).kern == max(stat(i).kern(:)));   % take slices through max
    
    resp_sf = stat(i).kern(ii(1),:);
    resp_ori = stat(i).kern(:,jj(1))';

    stat(i).resp_sf = resp_sf;
    stat(i).resp_ori = resp_ori;
    
    
    
    resp_sf = resp_sf - max(resp_sf)*.75; % clip tails...
    resp_sf(resp_sf<0) = 0;

    stat(i).ori_est = rad2deg(angle(sum(resp_ori.*exp(1i*ori*2*pi/180)))/2);
    if(stat(i).ori_est<0)
        stat(i).ori_est = stat(i).ori_est+180;
    end
    stat(i).sf_est  = 10^(sum(resp_sf.*log10(sf))/sum(resp_sf));
    
    stat(i).cv = 1-abs(sum(resp_ori.*exp(1i*2*pi*ori/180))./sum(resp_ori));
    
    hold on
    plot(log10(stat(i).sf_est),stat(i).ori_est,'k.','markersize',15);
    hold off;
    
    title(sprintf('Cell #%d',i));
        
 
    % try to reconstruct linear rf - decimate by 4

%     [xx,yy] = meshgrid((1:1920/4)-1920/8,(1:1080/4)-1080/8);
%     rf = zeros(size(xx));
%     sp = sper_pix/4;
%     for m=1:max_ori
%         for j=1:max_sper
%             for k=1:max_sphase
%                 if(N(m,k,j)>0)
%                     sth = sind(ori(m));
%                     cth = cosd(ori(m));
%                     stim = cos (2*pi * (cth*xx + sth*yy) / sp(j) + sphase(k));
%                     rf = rf + R(m,k,j,stat(i).tmax,i)*stim;
%                 end
%             end
%         end
%     end
% 
%     subplot(1,2,2)
%     imagesc(rf)
%     axis off;
    
   pause(0.05);
   %pause;
end

ori = 0:10:170;
sper = 1920./(1.31.^(0:12));   % spatial period in pixels
sper = sper / 9.6;             % 9.6 pixels per degree near center of screen
sf = 1./sper;               % cycles/deg
sphase = linspace(0,360,9);
ph = sphase(1:end-1);


save([fname '.orisf'],'r','stat','ori','sf','ph','-v7.3');
disp('Done!');



