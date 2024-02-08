function localizer(animal,unit,expt)

display_param;

open_ovserver
open_sbserver

nreps = 4; % (default is 4 for short version, do 14 for long versions)

rtime = 9*5*nreps;%12*7*4; %% run time in seconds

send_sbserver(sprintf('A%s',animal));
send_sbserver(sprintf('U%s',unit));
send_sbserver(sprintf('E%s',expt));

set_logfile(sprintf('%s_%s_%s.log',animal,unit,expt));

delete_object(-1);

add_object(0,17);            % Localizer...

set_trial(0);

send_sbserver(sprintf('G')); % tell microscope to start sampling
pause(10);                   % let it go for ... resonant mirror warm up

loop(rtime);  % run it for rtime min 
pause(rtime+18);

send_sbserver('S');     % stop

close_ovserver
close_sbserver

pause(2);
analyze_localizer(animal,unit,expt)
