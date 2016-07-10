function [H, delay] = Hgrpdelay(H, f)
%% Calculate and remove group delay from transfer function H
phiH = unwrap(angle(H));
df = abs(f(1)-f(2));
delay = -diff([phiH 0])/df; % in s
delay = interp1(f, delay, 0);

H = H.*exp(1j*2*pi*f*delay);