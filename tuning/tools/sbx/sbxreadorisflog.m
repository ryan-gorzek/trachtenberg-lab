function log = sbxreadorisflog(fname)

fn = [fname '.log_31'];
fid = fopen(fn,'r');
l = fgetl(fid);
k=0;

log = {};
while(l~=-1)
    k = str2num(l(2:end));
    log{k+1} = [];
    e = [];
    l = fgetl(fid);
    while(l~=-1 & l(1)~='T')
        e = [e ; str2num(l)];
        l = fgetl(fid);
    end
    log{k+1} = table(e(:,1),e(:,2),e(:,3),e(:,4),...
    'VariableNames',{'frame' 'ori' 'sphase' 'sper'});
end

sbx = load(fname);

% just in case TTL1 is not connected
% scanbox_frame = sbx.info.frame(2:end-1);

if isfield(sbx.info.evt, 'photodiode_off')
    scanbox_frame = sbx.info.evt.photodiode_off.frame;
elseif isfield(sbx.info.evt,'frame')
    scanbox_frame = sbx.info.evt.frame.frame;
elseif isfield(sbx.info.evt,'reward')
    scanbox_frame = sbx.info.evt.reward.frame;
else
    scanbox_frame = sbx.info.evt.sync.frame;
end

if numel(scanbox_frame) > size(log{1}, 1)
    scanbox_frame = scanbox_frame(1:end-1);
end

% % detect missing TTLs and fill in...
% 
% d = diff(scanbox_frame);
% du = unique(d);
% 
% while ~all(du<7)
%     idx = find(d>=7);
%     idx = idx(1);
%     scanbox_frame = [scanbox_frame(1:idx) ; round((scanbox_frame(idx)+scanbox_frame(idx+1))/2) ; scanbox_frame(idx+1:end)];
%     d = diff(scanbox_frame);
%     du = unique(d);
% end

% add sbxframe to all logs...

k = 1;
for j = 1:length(log)
    if ~isempty(log{j})
%         ov_frame = log{j}.frame;
%         fit = polyfit(ov_frame,scanbox_frame(k:k+length(ov_frame)-1),1);           % from ov to sbx
%         k = k+length(ov_frame);
%         t = table(floor(polyval(fit,ov_frame)),'VariableName',{'sbxframe'});
        t = table(scanbox_frame,'VariableName',{'sbxframe'});
        log{j} = [log{j} t];
    end
end

