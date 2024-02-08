
function r = battery2_ori_tf(animal,unit,expt,xpos,ypos,radius,sf)

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

add_object(0,6);            % creates a grating... and adds some default values
set_field('xpos',0,xpos,1); 
set_field('ypos',0,ypos,1);  
set_field('tper',0,20,1);   % 3Hz
set_field('sper',1,opt_sper,1);  % 0.04 cycles per deg
set_field('radius',0,round(radius*pixPerDeg),1);
set_field('contrast',1,0.9,1);


% create the matrix of orientations (directions) and spatial frequencies 

ori_points = 0:45:360; ori_points(end) = [];
tper_points = [8    12    20    30    48    60   120];
loop_points = [123, 123, 123, 123, 147, 123, 123]; % number of frames to display for each tper

[ori,tper]=meshgrid(ori_points,tper_points);

param = [ori(:) tper(:)]; % all the parameters!

send_sbserver(sprintf('G')); % tell microscope to start sampling
pause(10);                   % let it go for ... resonant mirror warm up

total_rpts = 8;             % total # of repeats 

perms = cell(1,total_rpts);

for(rpt=1:total_rpts)             % for each repeat

    p = randperm(size(param,1));  % new permutation 
    perms{rpt} = p;

    for(k = 1:length(p))

        anOri = param(p(k),1);
        aTper = param(p(k),2);
        
        set_field('th',1,anOri,1);
        set_field('tper',0,aTper,1);
               
        loop(loop_points(tper_points == aTper),1);         % present the stimulus for 1 sec... 
        pause(3.5);            % 1 sec blank
        
    end
end

pause(1);

send_sbserver('S');     % stop the microscope...

close_ovserver    % close communication channels...
close_sbserver

fn = [ animal '_' unit '_' expt '_p.mat'];

save(fn,'param','perms'); % save order of parameters

