
% real time analysis...

d = dir('*.sbx'); % get the file name
fname = d(1).name(1:end-4);

load([fname '.mat']); % contains the syncs and it is loaded into an info variable in the global space

frame_on  = info.evt.stim_on.frame;   % stim onsets
frame_off = info.evt.stim_off.frame;  % stim offsets

load([fname '_realtime.mat']);        % loads real time signals (raw fluorescence (rtdata frames x ncell matrix) and deconvolved spikes (stdata) )

% rt = bsxfun(@minus,mean(rtdata),rtdata);    % mean subtract
rt(:,isnan(mean(rtdata))) = [];             % remove bad ROIs

srt = medfilt1(rt,101);
rt = srt-rt;

ncells = size(rt,2);                  % how many cells
nstim = length(frame_on);             % total number of stim presented

vnorm = zeros(nstim,ncells);

for i = 1:nstim                       % for each stim index
    v = mean(rt(frame_on(i):frame_off(i),:));
    vnorm(i,:) = v;
end

load([fname '_p.mat']);              % parameters of stimuli (param) and order of presentation (perms)
                                     % perms is a cell array where each
                                     % entry is the permutation of stimuli
                                     % for that block ...
idx = horzcat(perms{:});

M = max(idx);                       % the maximum number of stimuli 
navg = zeros(M,ncells);             % and the average (across repeats) response to each 

for k = 1:M
    j = find(idx==k);               % find all the times we presented stim k
    navg(k,:) = mean(vnorm(j,:));   % average the response
end

% this is reconstructing the stim matrix because I didn't save the matrices
% in the _p file.

p1 = unique(param(:,1)); %#ok<*USENS> 
p2 = unique(param(:,2));

nparam1 = length(p1);
nparam2 = length(p2);
[pp1,pp2] = meshgrid(p1,p2);

navg = reshape(navg,[size(pp1) ncells]);     % reshape so it has the right matrix form

sfpop = mean(navg,2);               % average across all orientations to generate a mean spatial frequency tuning curve        

plot(log10(p2),mean(navg,2),'-o')   % plot it...
xlabel('Log10(spatial freq)');
ylabel('Population response');

p=polyfit(log10(p2),mean(navg,2),2);    % fit a parabola (poly deg 2)
opt_sf = 10^(-p(2)/(2*p(1)));           % -b/(2a)
title(['Optimal sfreq = ' num2str(opt_sf,2)]);  % choose that as the center (optimal) sf
figure(gcf)

% % lets look at the temporal response for the spikes
% 
% savg = zeros(1,45);
% for i=1:nstim
%     savg = savg + vecnorm(rt(frame_on(i)-15 : frame_on(i)+29, :),2,2)';
% end

