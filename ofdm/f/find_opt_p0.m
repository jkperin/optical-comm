% Find initial point of the quantizer to maximize resolution
% Given the input signal (x), resoltution (delta), and range (Xtot), it
% calculates the initial point of the quantizer so that the most part of
% the signal are inside the quantization range
function n0 = find_opt_p0(x, delta, Xtot)

Nbins = 100;

[fx, xl] = hist(x, Nbins);
fx = fx/trapz(fx); 

n0 = min(ceil(min(x)/delta), 0);
x0 = n0*delta;
S = trapz(fx(xl >= x0  & xl <= x0 + Xtot));
Sprev = 0;
while Sprev < S
    Sprev = S;
    n0 = n0 + 1;
    x0 = n0*delta;
    S = trapz(fx(xl >= x0  & xl <= x0 + Xtot));
end

n0 = n0 - 1;

