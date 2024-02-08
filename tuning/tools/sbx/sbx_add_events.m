function r = sbx_add_events(info)

evt = unique(info.event_id);
if ~isfield(info,'evt_names')
    r = [];
end

evt_names = info.evt_names;

t = (info.frame*info.sz(1)+info.line)/info.resfreq;    % in seconds

N = length(evt);
for i = 1:N
    idx = info.event_id==evt(i);
    n = evt_names{round(log2(double(evt(i))))+1};
    cmd = sprintf('r.%s.time = t(idx);',n); eval(cmd);
    cmd = sprintf('r.%s.frame = info.frame(idx);',n); eval(cmd);
    cmd = sprintf('r.%s.line  = info.line(idx);',n); eval(cmd);
end
