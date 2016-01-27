% scprop.m   Joseph M. Kahn  9-5-11
% Models propagation of optical field through a scalar system (such as a WSS or SMF with negligible PMD and PDL).
function Efilt = scprop(Ein,Hfilt)
Efilt = [ifft(fft(Ein(1,:)).*Hfilt); ifft(fft(Ein(2,:)).*Hfilt)];
end





