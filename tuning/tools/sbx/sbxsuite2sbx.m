function sbxsuite2sbx(smat,fn)

% Convert a suite2p segmentation to Scanbox format

% z = sbxread(fn,0,1);
% z = squeeze(z(1,:,:));
% sz = size(z);

sz = [512 796];
mask = zeros(sz);

s2p=load(smat);

k = 1;
spks = [];
sig = [];
np = [];
m = [];
n = [];

for i=1:length(s2p.stat)
    if(s2p.iscell(i,1))
        mask(sub2ind(sz,s2p.stat{i}.ypix,s2p.stat{i}.xpix)) = k;
        spks(:,k) = s2p.spks(i,:)';
        sig(:,k)  = s2p.F(i,:);
        np(:,k)   = s2p.Fneu(i,:)';
        k = k+1;
    end
end
    
if s2p.ops.nchannels == 1
    m = s2p.ops.meanImg;
    p = s2p.ops.max_proj;
    
elseif s2p.ops.nchannels == 2
    m = s2p.ops.meanImg;
    p = s2p.ops.max_proj;
    n = s2p.ops.meanImg_chan2;
    r0 = s2p.ops.meanImg_chan2_corrected;
    a = s2p.ops.yrange;
    b = s2p.ops.xrange;
    r = r0((a(1)+1):a(2),(b(1)+1):b(2));
    M = (r-min(r(:)))/(max(r(:))-min(r(:)));
    x = adapthisteq(M);
    x = single(x);
    x = (x-min(x(:)))/(max(x(:))-min(x(:)));
    merge = imfuse(x,p,'falsecolor','Scaling','independent','ColorChannels',[1 2 0]);
end
if ~isempty(n)
    save([fn '.align'],'m','n','p','r','x','merge');
    save([fn '.segment'],'mask');
    save([fn '.signals1'],'spks','sig','np');
else
    save([fn '.align'],'m','p');
    save([fn '.segment'],'mask');
    save([fn '.signals1'],'spks','sig','np');
end


