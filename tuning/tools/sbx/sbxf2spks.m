function sbxf2spks(fname)

load([fname '.signals1'],'-mat'); % load signals

ncell = size(sig,2);
spks = zeros(size(sig,1),ncell);

for i = 1:ncell
    
    spks(:,i) = deconv(sig(:,i), [1.5    1.8503    0.2958    9.3894]);
    
end
                
    save([fname '.signals'],'sig','np','spks');     
 




function z = deconv(y,x)

s = x(1);   % sigma
th = x(2);  % theta
b = x(3);   % beta
a = x(4);   % alpha

nsamp = size(y,1);

% Odd filter

t = 0:nsamp-1;
w = t.*exp(-t.^2 / (2*s^2));
w(2:end) = (w(2:end)-w(end:-1:2));
w = -w';
w = w/norm(w);

% Even filter

w0 = zeros(nsamp,1);
w0 = exp(-t.^2 / (2*s^2));
w0(2:end) = (w0(2:end)+w0(end:-1:2));
w0 = w0';
w0 = w0/norm(w0);
 
% Filtered signals

wf0 = fft(w0);
xf0 = real(ifft(fft(y(1:nsamp)).*wf0));
xf0 = zscore(xf0);

wf = fft(w);
xf = real(ifft(fft(y(1:nsamp)).*wf));
xf = zscore(xf);

% Of course one can combine the filters first and convolve once...
% but for historical reasons I kept them separate.

% Linear combination of filtered signals

xf = cosd(a)*xf+sind(a)*xf0;

% Output nonlinearity

z(1:nsamp) = (xf-th).^b .* (xf>=th);


