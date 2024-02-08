
function r = battery1_ori_sf(animal,unit,expt,xpos,ypos,radius)

% xpos 
% ypos from localizer... (500,1000)
% radius - size of window in DEG!!!

clear param;

display_param;

open_ovserver   % open connection to the stimulus computer
open_sbserver   % open connnection to the microscope 

send_sbserver(sprintf('A%s',animal)); % tell the microscope the animal name...
send_sbserver(sprintf('U%s',unit));   % the ROI name...
send_sbserver(sprintf('E%s',expt));   % the experient number

% setup  grating

delete_object(-1);

add_object(0,6);            % creates a grating... and adds some default values
set_field('xpos',0,xpos,1); 
set_field('ypos',0,ypos,1);  
set_field('tper',0,60,1);   % 3Hz
set_field('sper',1,round(pixPerDeg/.04),1);  % 0.04 cycles per deg
set_field('radius',0,round(radius*pixPerDeg),1); 
set_field('contrast',1,0.9,1);


% create the matrix of orientations (directions) and spatial frequencies 



ori_points = 0:45:360; ori_points(end) = [];
sf_points = logspace(log10(0.01),log10(0.30),7);

[ori,sf]=meshgrid(ori_points,sf_points);

param = [ori(:) sf(:)]; % all the parameters!

send_sbserver(sprintf('G')); % tell microscope to start sampling
pause(10);                   % let it go for ... resonant mirror warm up

total_rpts = 8;             % total # of repeats 

perms = cell(1,total_rpts);

for(rpt=1:total_rpts)             % for each repeat

    p = randperm(size(param,1));  % new permutation 
    perms{rpt} = p;

    for(k = 1:length(p))

        anOri = param(p(k),1);
        aSpatialPer = pixPerDeg / param(p(k),2);
        
        set_field('th',1,anOri,1);
        set_field('sper',1,aSpatialPer,1);
               
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

