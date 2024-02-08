
d = dir('*.sbx'); % get the file name
fname = d(1).name(1:end-4);

load([fname '.mat']); % contains the syncs and it is loaded into an info variable in the global space
load([fname '_quadrature.mat']); % contains the wheel data
quad_diff = abs(diff(quad_data));

figure; tiledlayout(2, 1);
nexttile; plot(quad_data, "LineWidth", 2);
nexttile; plot(quad_diff, "LineWidth", 2);

frame_on  = info.evt.stim_on.frame;   % stim onsets
frame_off = info.evt.stim_off.frame;  % stim offsets

nstim = numel(frame_on);

load([fname '_p.mat']);              % parameters of stimuli (param) and order of presentation (perms)
                                     % perms is a cell array where each
                                     % entry is the permutation of stimuli
                                     % for that block ...
idx = horzcat(perms{:});
M = max(idx);

stim_quad = zeros(nstim, 1);
for s = 1:nstim
    stim_quad(s) = mean(quad_diff(frame_on(s) : frame_off(s)));
end

navg_quad = zeros(M, 1);
for k = 1:M
    j = find(idx==k);               % find all the times we presented stim k
    navg_quad(k) = mean(stim_quad(j));
end

p1 = unique(param(:,1)); %#ok<*USENS> 
p2 = unique(param(:,2));
nparam1 = numel(p1);
nparam2 = numel(p2);
[pp1,pp2] = meshgrid(p1, p2);
navg_quad = reshape(navg_quad, size(pp1));

figure; tiledlayout(1, 3);
nexttile; imagesc(navg_quad); colorbar;
xticks(1:nparam1); xticklabels(round(p1, 3)); xlabel("Direction");
yticks(1:nparam2); yticklabels(round(p2, 3)); ylabel("Spatial Freq.")
nexttile; plot(round(p1, 2), mean(navg_quad, 1), "-o", "LineWidth", 2);
xticks(round(p1, 3)); xlabel("Direction"); ylabel("abs(\DeltaQuadrature)");
axis tight;
nexttile; plot(round(p2, 2), mean(navg_quad, 2), "-o", "LineWidth", 2);
xticks(round(p2, 3)); xlabel("Spatial Freq."); ylabel("abs(\DeltaQuadrature)");
axis tight;
