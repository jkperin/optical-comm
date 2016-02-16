% freqshift.m   Joseph M. Kahn      12/6/11
% Shifts a vector signal by frequency fshift.
% x: input signal (two-row matrix)
% t: time (one-row vector)
function xshift = freqshift(x,t,fshift)
%xshift = [x(1,:).*exp(sqrt(-1)*2*pi*fshift*t); x(2,:).*exp(sqrt(-1)*2*pi*fshift*t)];
xshift = x.*repmat(exp(1j*2*pi*fshift*t),size(x,1),1);
end
