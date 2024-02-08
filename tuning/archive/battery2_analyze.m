
d = dir('*.sbx');
fname = d(1).name(1:end-4);

load([fname '.mat']);

frame_on = info.evt.stim_on.frame; % stim onsets
frame_off = info.evt.stim_off.frame; % stim offsets

load([fname '_realtime.mat']);

rt = bsxfun(@minus,mean(rtdata),rtdata);
rt(:,isnan(mean(rtdata))) = []; % remove bad ROIs

nstim = length(frame_on);

vnorm = zeros(1,nstim);

for i = 1:nstim
    v = rt(frame_on(i):frame_off(i),:);
    v0 = rt(frame_on(i)-1);
    v = bsxfun(@minus,v,v0);
    vnorm(i) = mean(vecnorm(v,2,2));
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

tfpop = mean(navg,2);
tfreq = 60./p2;

plot(tfreq,mean(navg,2),'-o')
xlabel('Temporal freq (Hz)');
ylabel('Population response');

% p=polyfit(log10(p2),mean(navg,2),2);
% opt_sf = 10^(-p(2)/(2*p(1)));
% title(['Optimal tf = ' num2str(opt_sf,2)]);


