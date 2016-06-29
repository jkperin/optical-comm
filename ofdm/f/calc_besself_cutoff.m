%% Function used to determine the conversion factor between cutoff frequency and w0, which is the input parameter of the matlab function besself
% The wo parameter in Bessel filter design is the frequency up to which the
% filter's group delay is approximately constant. 

% Example:
% [x, fval] = fzero(@(x) calc_besself_cutoff(5, 1, x) - 0.5, 1.65)
% x =
%    1.621597623980423
% 
% fval =
%      3.330669073875470e-16
% Calculates wc2w0 (i.e., x) that makes f3dB = 1 the cutoff frequency of
% a 5th-order Bessel filter (i.e., |h(f3dB)|^2 == 0.5)
% This result is independent of the actual frequency f3dB.

function h = calc_besself_cutoff(order, f3dB, wc2w0)

[b,a] = besself(order, 2*pi*f3dB*abs(wc2w0));

h = polyval(b, 1j*2*pi*f3dB)/polyval(a, 1j*2*pi*f3dB);

h = abs(h).^2;
