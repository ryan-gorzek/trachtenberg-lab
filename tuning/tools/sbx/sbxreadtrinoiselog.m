function log = sbxreadtrinoiselog(fname,mform)

global info

fn = [fname '.log_59'];
fid = fopen(fn,'r');
l = fgetl(fid); 
log = {};

load(fname); % info

r = sbx_add_events(info);

try
  f = r.frame.frame;    % get the frames
catch
  f = r.sync.frame;  
end

m = 1;
while(l~=-1)
    k = str2num(l);
    t = k(1);
    k = reshape(k(2:end),3,[])';
    stim = full(sparse(k(:,1)+1,k(:,2)+1,k(:,3)));
    log{end+1} = {t stim f(m)};
    l = fgetl(fid);
    m = m+1;
end




