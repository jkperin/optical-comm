%% Frequency locked loop for QPSK

clear, clc, close all

addpath ../../f/

% Parameters
Nsymb = 2^10;
Mct = 15;
N = Nsymb*Mct;
Rs = 56e9;
fs = Rs*Mct;
df = fs/N;
f = -fs/2:df:fs/2-df;
dt = 1/fs;
t = 0:dt:(N-1)*dt;
ts = t;
Laser = laser(1310e-9, 0, -150, 200e3);

% Simulation control
phase_noise = ~true;
awgn_noise = ~true;
frequency_offset = true;

% Generate symbols
constellation = pi/4*[1 3 5 7];
data = randi([1 4], [Nsymb 1]);
x = qammod(data-1, 4, 0, 'Gray');
xi = real(x);
xq = imag(x);
di = constellation(data).';
dq = constellation(data).';

% Pulse shape
xi = reshape(repmat(xi, Mct, 1), [], N);
di = reshape(repmat(di, Mct, 1), [], N);
xq = reshape(repmat(xq, Mct, 1), [], N);
dq = reshape(repmat(dq, Mct, 1), [], N);
x = xi + 1j*xq;

% include phase noise
phin = 0;
if phase_noise
    [~, phin] = Laser.addPhaseNosie(x, fs);  
end

% include additive noise
varN = 0.0001;
nn = 0;
if awgn_noise
    nn = sqrt(varN/2)*randn(size(x)) + 1j*sqrt(varN/2)*randn(size(x));
end

% include frequency offset
foff = 3e9;
foff = double(frequency_offset)*foff
xi = cos(2*pi*foff*t + di + phin);
xq = sin(2*pi*foff*t + dq + phin);

xn = xi + 1j*xq + nn;
Hfilt = design_filter('Bessel', 5, 20e9/(fs/2));
% 
% xi = real(ifft(fft(real(xn)).*ifftshift(Hfilt.H(f/fs))));
% xq = real(ifft(fft(imag(xn)).*ifftshift(Hfilt.H(f/fs))));
% xn = ifft(fft(xn).*ifftshift(Hfilt.H(f/fs)));
% xn = xi + 1j*xq;

% Loop filter
csi = 1/sqrt(2);                                                    % damping coefficient of second-order loop filter
wn = 2*pi*1e9;                                                    % relaxation frequency of second-order loop filter: optimized using optimize_PLL.m
Kdc = 1;
nums = Kdc*[2*csi*wn wn^2];
dens = [1 0 0]; % descending powers of s
% %
% G = tf(nums, dens);
% H = G/(1 + G);
% [nums, dens] = tfdata(H);
% nums = cell2mat(nums);
% dens = cell2mat(dens);

[numz, denz] = impinvar(nums, dens, fs);
numLen = length(numz);
denLen = length(denz);
% % 
% y = xqd.*xin - xid.*xqn;
% yf = real(ifft(fft(y).*fft(fBW*(sinc(t*fBW).^2))));

xr = zeros(size(x));
y = zeros(size(x));
yf = zeros(size(x));
ydiff = zeros(size(x));
D = 0;
for t = max([numLen, D, length(Hfilt.h) 50])+1:length(x)
    xr(t) = exp(1j*(-yf(t-1)))*xn(t);
    xr(t) = xn(t);
    
    xi(t) = real(xr(t));
    xq(t) = imag(xr(t));
    xid(t) = sign(xi(t));
    xqd(t) = sign(xq(t));
        
    % Costas loop
    y(t) = xid(t-D)*xq(t) - xqd(t-D)*xi(t);
%     y(t) = -xi(t)*xq(t-D) + xq(t)*xi(t-D);

    % XOR
    comp = (abs(xi(t)) < abs(xq(t)));
    tmp = not(xor(xi(t) >= 0, xq(t) >= 0));
    y(t) = not(xor(tmp, comp));
    y(t) = 2*y(t) - 1;

    yf1(t) = sum(Hfilt.h.*y(t:-1:t-length(Hfilt.h)+1));
    yf1(t) = y(t);
    
    ydiff(t) = abs(yf1(t) - yf1(t-1));
        
    yf(t) = sum(numz.*ydiff(t:-1:t-numLen+1)) - sum(yf(t-1:-1:t-denLen+1).*denz(2:end));
end

figure
subplot(211), hold on
plot(y)
% plot(yf1)
% plot(ydiff)
subplot(212), hold on
plot(yf)
plot(2*pi*foff*ts)


scatterplot(xn)
title('input')
scatterplot(xr(end-1000:end))
title('output')

% H = tf(numz, denz, 1/fs);
% woff = 2*pi*foff;
% t = 0:1/fs:2^10/fs;
% phin = pi/5;
% 
% u = sin(woff*t + phin);
% yy = lsim(H, u, t);
% 
% ytheory = (wn/woff)^2*(-sin(woff*t+phin) + sin(phin)+woff*t*cos(phin))...
%     -2*csi*wn/woff*(cos(woff*t+phin)-cos(phin));
% 
% figure
% plot(t, u, t, yy, t, ytheory, '--')

% wvec = woff;
% for k = 2:20;
%     wo = wvec(k-1) - 
%     wvec(k) = wvec(k-1) + wn^2/wvec(k-1);
% end
% 
% figure, plot(wvec/(2*pi))
    
    
    


