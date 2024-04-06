close all;
clearvars;
clc;

load('data_all'); % load pre-trained solutions

%%
nsol = length(data_all.err); % number of pre-trained solutions
nc = 4; % number of cell-types in the recurrent sub-circuit (exclusing feedforward input)
nb = 4; % number of sub-circuits
nc2 = nc*nb; % number of units across recurrent sub-circuits

is0 = 2; % select center stimulus
doo = 0.05; % optogenetic strength

new_max = 1.; % max activity allowed (otherwise classified as unstable)


%% create array with contrast values
nco = 11; 
contrast_c = linspace(0.,1.,nco);

%% compute responses to a center stimulus in absence of optogenetic perturbations (control)


fixed_points = zeros(nsol,nco,nc2);
num_iters = zeros(nsol, nco);
dynpars = data_all.dynpars;
for ico = 1:nco
    fprintf('%d',ico)
    cc = contrast_c(ico);
    for isol = 1:nsol
        if mod(isol,10) == 0
            fprintf('.')
        end
        rt = nan(nc2,1);
        rt(:,1) = squeeze(data_all.X(isol,1:nc2,is0)); % starting point (for a center stimulus)
        h = cc*squeeze(data_all.X(isol,nc2+1:end,is0))'; % L4 input (for a center stimulus)
        W = squeeze(data_all.W(isol,:,:)); % W matrix
        dr = 1;
        it = 0;
        while it < dynpars.nt && dr > dynpars.eps_dr && max(rt(:)) < new_max % stop if timeout | converged | unstable (high activity)
            it = it + 1;
            r0 = rt(:,it);
            r1 = r0 + (-r0 + subplus(W*[r0; h]).^2 ).*dynpars.dt;
            rt = cat(2, rt, r1);
            dr = sum((r1(:) - r0(:)).^2);

        end
        fixed_points(isol,ico,:) = r1;
        num_iters(isol, ico) = it;
    end
    fprintf('\n')
end


%% optogenetic perturbation

r_opto_pv_t = cell(nco,nsol);
r_opto_pv = nan(nco,nc2,nsol);
time_steps_pv = nan(nco,nsol);

r_opto_vip_t = cell(nco,nsol);
r_opto_vip = nan(nco,nc2,nsol);
time_steps_vip = nan(nco,nsol);


for ico = 1:nco
    fprintf('%d',ico)
    cc = contrast_c(ico);
    for iu = 1:2
        if iu == 1 % PV excitation
            iopto = 2;
            
        else
            iopto = 4; % VIP excitation
            
        end
        
        Wo = zeros(nc,1);
        Wo(iopto) = 1;
        Wo = repmat(Wo,nb,1);
        
        Oh = Wo.*doo; % optogenetic input
        
        dynpars = data_all.dynpars;
        for isol = 1:nsol
            if mod(isol,10) == 0
                fprintf('.')
            end
            rt = nan(nc2,1);
            rt(:,1) = squeeze(data_all.X(isol,1:nc2,is0)); % starting point (for a center stimulus)
            h = cc*squeeze(data_all.X(isol,nc2+1:end,is0))'; % L4 input (sor a center stimulus)
            W = squeeze(data_all.W(isol,:,:)); % W matrix
            dr = 1;
            it = 0;
            while it < dynpars.nt && dr > dynpars.eps_dr && max(rt(:)) < new_max
                it = it + 1;
                r0 = rt(:,it);
                r1 = r0 + (-r0 + subplus(W*[r0; h] + Oh).^2 ).*dynpars.dt;
                rt = cat(2, rt, r1);
                dr = sum((r1(:) - r0(:)).^2);
            end
            
            if iu == 1
                if it < dynpars.nt && max(rt(:)) < dynpars.max_x
                    r_opto_pv(ico,:,isol) = r1;
                end
                r_opto_pv_t{ico,isol} = rt;
                time_steps_pv(ico,isol) = it;
            else
                if it < dynpars.nt && max(rt(:)) < dynpars.max_x
                    r_opto_vip(ico,:,isol) = r1;
                elseif max(rt(:)) < dynpars.max_x
                    disp("Unstable...");
                end
                r_opto_vip_t{ico,isol} = rt;
                time_steps_vip(ico,isol) = it;
            end
            
        end
        fprintf('\n')
    end
end

sol_exist = find(squeeze(~isnan(sum(r_opto_vip(:,1,:),1)) & ~isnan(sum(r_opto_pv(:,1,:),1)))); % find stable and existing solutions

%% average activity across solutions
exc_m = squeeze(mean(fixed_points(sol_exist,:,1),1)); % average excitatory activity across solutions (control)
pv_m = squeeze(mean(r_opto_pv(:,1,sol_exist),3)); % average excitatory activity across solutions (PV excitation)
vip_m = squeeze(mean(r_opto_vip(:,1,sol_exist),3)); % average excitatory activity across solutions (VIP excitation)
%%
figure;
hold on
plot(contrast_c,exc_m,'k-','Linewidth',1.5)
plot(contrast_c,pv_m,'b-','Linewidth',1.5)
plot(contrast_c,vip_m,'g-','Linewidth',1.5)
legend({'Control','PV activation','VIP activation'},'fontsize',18,'Location','NorthWest')
set(gca,'fontsize',18)
xlabel('Contrast','fontsize',18)
ylabel('Response','fontsize',18)
axis square


figure;
hold on
plot(log(contrast_c),log(exc_m),'k-','Linewidth',1.5)
plot(log(contrast_c),log(pv_m),'b-','Linewidth',1.5)
plot(log(contrast_c),log(vip_m),'g-','Linewidth',1.5)
legend({'Control','PV activation','VIP activation'},'fontsize',18,'Location','NorthWest')
set(gca,'fontsize',18)
xlabel('Log(contrast)','fontsize',18)
ylabel('Log(response)','fontsize',18)
axis square

