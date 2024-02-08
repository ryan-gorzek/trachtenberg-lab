
d = dir('*.sbx');
fname = d(1).name(1:end-4);

load([fname '.mat']);

frame_on = info.evt.stim_on.frame; % stim onsets
frame_off = info.evt.stim_off.frame; % stim offsets

load([fname '_realtime.mat']);

% rt = bsxfun(@minus,mean(rtdata),rtdata);
rt = stdata;
rt(:,isnan(mean(rtdata))) = []; % remove bad ROIs

nstim = length(frame_on);

vnorm = zeros(1,nstim);

for i = 1:nstim
    v = rt(frame_on(i):frame_off(i),:);
    v0 = rt( (frame_on(i)-3) : (frame_on(i)-1) , :);
    vsub = bsxfun(@minus,v,mean(v0,1));
    vnorm(i) = mean(vecnorm(vsub,2,2));
end

load([fname '_p.mat']);

idx = horzcat(perms{:});

M = max(idx);

navg = zeros(1,M);

for k = 1:M
    j = find(idx==k);
    navg(k) = mean(vnorm(j));
end

p1 = unique(param(:,1)); %#ok<*USENS> 
p2 = unique(param(:,2));
nparam1 = length(p1);
nparam2 = length(p2);
[pp1,pp2] = meshgrid(p1,p2);

navg = reshape(navg,size(pp1));

sfpop = mean(navg,2);

plot(log10(p2),mean(navg,2),'-o')
xlabel('Log10(spatial freq)');
ylabel('Population response');

p=polyfit(log10(p2),mean(navg,2),2);
opt_sf = 10^(-p(2)/(2*p(1)));
title(['Optimal sfreq = ' num2str(opt_sf,2)]);


