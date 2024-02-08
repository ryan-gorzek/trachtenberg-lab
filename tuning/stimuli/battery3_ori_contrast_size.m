
function r = battery3_ori_contrast_size(animal,unit,expt,xpos,ypos,sf,tf)

clear param;

display_param;

open_ovserver   % open connection to the stimulus computer
open_sbserver   % open connnection to the microscope 

send_sbserver(sprintf('A%s',animal)); % tell the microscope the animal name...
send_sbserver(sprintf('U%s',unit));   % the ROI name...
send_sbserver(sprintf('E%s',expt));   % the experient number

% setup  grating

delete_object(-1);

opt_sper = round(pixPerDeg/sf);
opt_tper = 60/tf;

add_object(0,6);            % creates a grating... and adds some default values
set_field('xpos',0,xpos,1); 
set_field('ypos',0,ypos,1);  
set_field('tper',0,opt_tper,1);  
set_field('sper',1,opt_sper,1);  % 0.04 cycles per deg
set_field('radius',0,round(10*pixPerDeg),1); 
set_field('contrast',1,0.9,1);

% create the matrix of orientations (directions) and spatial frequencies 

ori_points = 0:45:360; ori_points(end) = [];  % orientation (direction)
% c_points = [0  0.0500    0.0824    0.1357    0.2236    0.3684    0.6070    1.0000]; % contrast
c_points = [0 0.0794 0.1318 0.2188 0.3631 0.6026, 1.0000];
r_points = [2.5 5 10 15 20 25 30]; % radius in deg

[ori,cont,radius]=meshgrid(ori_points,c_points,r_points);

param = [ori(:) cont(:) radius(:)]; % all the parameters!

send_sbserver(sprintf('G')); % tell microscope to start sampling
pause(10);                   % let it go for ... resonant mirror warm up

total_rpts = 5;             % total # of repeats 

perms = cell(1,total_rpts);

for(rpt=1:total_rpts)             % for each repeat

    p = randperm(size(param,1));  % new permutation 
    perms{rpt} = p;

    for(k = 1:length(p))

        anOri = param(p(k),1);
        aCont = param(p(k),2);
        aRad =  round(param(p(k),3)*pixPerDeg);
        
        set_field('th',1,anOri,1);
        set_field('contrast',1,aCont,1);
        set_field('radius',0,aRad,1);

        loop(63,1);             % present the stimulus for 1 sec... 
        pause(2);            % 1 sec blank
        
    end
end

pause(1);

send_sbserver('S');     % stop the microscope...

close_ovserver    % close communication channels...
close_sbserver

fn = [ animal '_' unit '_' expt '_p.mat'];

save(fn,'param','perms'); % save order of parameters

