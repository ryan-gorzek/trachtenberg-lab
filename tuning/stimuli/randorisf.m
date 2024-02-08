
function randorisf(animal,unit,expt,xpos,ypos,radius,nori,nsper,cont,length)

display_param;

open_ovserver
open_sbserver

send_sbserver(sprintf('A%s',animal));
send_sbserver(sprintf('U%s',unit));
send_sbserver(sprintf('E%s',expt));

set_logfile(sprintf('%s_%s_%s.log',animal,unit,expt));

% setup  random grating object

delete_object(-1);
add_object(0,31);            % Hartley


set_field('tper',0,15,1);   % 4Hz; the monitor refresh rate is 60hz, so 60/15=4;
set_field('contrast',1,cont,1); %acceptable values run from 0 to 1
set_field('max_ori',0,nori,1); %number of orientations
set_field('max_sper',0,nsper,1); %displayWidth./ncycles); spatial period
set_field('radius',0,radius,1); %in pixels
set_field('xpos',0,xpos,1); %center of screen is 960
set_field('ypos',0,ypos,1); %center of screen is 540
set_trial(0);

% % flashing...
set_field('ton',0,5,1);
set_field('fmode',0,0,1);

send_sbserver(sprintf('G')); % tell microscope to start sampling
pause(10);                   % let it go for ... resonant mirror warm up



loop(60*length);  % run it for length min 
pause(60*length+8);    % and another 14 sec...

send_sbserver('S');     % stop

close_ovserver
close_sbserver


