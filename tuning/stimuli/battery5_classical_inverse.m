
function r = battery5_classical_inverse(animal,unit,expt,xpos,ypos,sf,tf)

clear param;

display_param;

open_ovserver   % open connection to the stimulus computer
% open_sbserver   % open connnection to the microscope 

% send_sbserver(sprintf('A%s',animal)); % tell the microscope the animal name...
% send_sbserver(sprintf('U%s',unit));   % the ROI name...
% send_sbserver(sprintf('E%s',expt));   % the experient number

% setup  grating

delete_object(-1);
set_ncomp(2);
add_object(0,39);            % creates a multigrating with 2 components

opt_sper = round(pixPerDeg/sf);
opt_tper = round(60/tf);
% opt_rad = round(pixPerDeg*rad);

set_array('sper',1,opt_sper,1,0);
set_array('th',1,0,1,0);
set_array('contrast',1,0.4,1,0);
% set_array('wt',0,2,1,0);
set_array('tper',0,opt_tper,1,0);
set_array('xpos',0,xpos,1,0);
set_array('ypos',0,ypos,1,0);
% set_array('radius',1,opt_rad,1,0);

set_array('sper',1,opt_sper,1,1);
set_array('th',1,0,1,1);
set_array('contrast',1,0.4,1,1);
% set_array('wt',0,4,1,1);
set_array('tper',0,opt_tper,1,1);
set_array('xpos',0,xpos,1,1);
set_array('ypos',0,ypos,1,1);
% set_array('radius',1,opt_rad,1,1);

% create the matrix of orientations (directions) and spatial frequencies 

ori_points = 0:90:360; ori_points(end) = [];  % orientation (direction)
c_points = [0 0.4]; % contrast (multiply by 2...)
r_points = 5:10:45;

[ori1,ori2,cont1,cont2,radius]=ndgrid(ori_points,ori_points,c_points,c_points,r_points);
param = [ori1(:) ori2(:) cont1(:) cont2(:) radius(:)]; % all the parameters!

% send_sbserver(sprintf('G')); % tell microscope to start sampling
pause(10);                   % let it go for ... resonant mirror warm up

total_rpts = 6;             % total # of repeats 

perms = cell(1,total_rpts);

for(rpt=1:total_rpts)             % for each repeat

    p = randperm(size(param,1));  % new permutation 
    perms{rpt} = p;

    for(k = 1:length(p))

        o1 = param(p(k),1);
        o2 = param(p(k),2);
        c1 = param(p(k),3);
        c2 = param(p(k),4);
        r = round(param(p(k),5)*pixPerDeg);

        set_array('th',1,deg2rad(o1),1,0);
        set_array('contrast',1,c1,1,0);
        set_array('radius',1,r,1,0);
        set_array('th',1,deg2rad(o2),1,1);
        set_array('contrast',1,c2,1,1);
        set_array('radius',1,r,1,1);
        
        loop(1);             % present the stimulus for 1 sec... 
        pause(2);            % 1 sec blank
        
    end
end

pause(1);

% send_sbserver('S');     % stop the microscope...

close_ovserver    % close communication channels...
% close_sbserver

fn = [ animal '_' unit '_' expt '_p.mat'];

save(fn,'param','perms'); % save order of parameters
