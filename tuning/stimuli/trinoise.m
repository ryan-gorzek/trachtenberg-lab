function trinoise(animal,unit,expt,xsize,ysize,len)

display_param;

open_ovserver
open_sbserver

send_sbserver(sprintf('A%s',animal));
send_sbserver(sprintf('U%s',unit));
send_sbserver(sprintf('E%s',expt));

set_logfile(sprintf('%s_%s_%s.log',animal,unit,expt));

% setup  random grating object

delete_object(-1);
add_object(0,59);            % trinoise

set_field('period',0,60,1);   % ISI
set_field('ton',0,15,1);      % duration of flash 
set_field('xsize',0,xsize,1); % checker size
set_field('ysize',0,ysize,1); 
set_field('maxrand',0,40,1);

dnforce_background;           % do not force background

set_trial(0);

send_sbserver(sprintf('G')); % tell microscope to start sampling
pause(10);                   % let it go for ... resonant mirror warm up

loop(60*len);       % run it for length min 
pause(60*len+16);    % and another 14 sec...

send_sbserver('S');     % stop

force_background;       % force background again

vars = get_fields;
save(sprintf('%s_%s_%s.vars',animal,unit,expt),'vars');

close_ovserver
close_sbserver

