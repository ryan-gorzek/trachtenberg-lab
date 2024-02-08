function [p, x, pv, opt] = fitcurve(x, y)

p = polyfit(x, y, 2);         % fit a parabola (poly deg 2)
opt = 10^(-p(2)/(2*p(1)));    % -b/(2a)
pv = polyval(p, x);
end