%% PLL for QPSK using XOR port for phase estimation

clear, clc, close all

addpath ../../f/

Nsymb = 2^10;
Mct = 15;
N = Nsymb*Mct;
Rs = 56e9;
fs = Rs*Mct;

Laser = laser(1310e-9, 0, -150, 2000e3);

phase_noise = true;
awgn_noise = true;
frequency_offset = true;

constellation = pi/4*[1 3 5 7];
data = randi([1 4], [Nsymb 1]);

% xi = cos(constellation(data));
% xq = sin(constellation(data));

x = qammod(data-1, 4, 0, 'Gray');
xi = real(x);
xq = imag(x);

xi = reshape(repmat(xi, Mct, 1), [], N);
xq = reshape(repmat(xq, Mct, 1), [], N);

x = xi + 1j*xq;

if phase_noise
    [xn, phin] = Laser.addPhaseNosie(x, fs);
else
    xn = x;
end

varN = 0.1;
if awgn_noise
    xn = xn + sqrt(varN/2)*randn(size(xn)) + 1j*sqrt(varN/2)*randn(size(xn));
end

foff = 10e9;
if frequency_offset   
    xn = freqshift(xn, 0:1/fs:(length(xn)-1)/fs, foff);
end

% Loop filter
csi = 1/sqrt(2);                                                    % damping coefficient of second-order loop filter
wn = 2*pi*1e9;                                                    % relaxation frequency of second-order loop filter: optimized using optimize_PLL.m
Kdc = 10;
nums = Kdc*[2*csi*wn wn^2];
dens = [1 0 0]; % descending powers of s

[numz, denz] = impinvar(nums, dens, fs);
numLen = length(numz);
denLen = length(denz);

xr = zeros(size(x));
y = zeros(size(x));
yf = zeros(size(x));
for t = Mct+numLen+1:length(x)
    xr(t) = exp(1j*(yf(t-4)))*xn(t);
    
    xi(t) = real(xr(t));
    xq(t) = imag(xr(t));
    
    xid(t) = 2*(xi(t) > 0) - 1;
    xqd(t) = 2*(xq(t) > 0) - 1;
    comp = 2*(abs(xi(t)) < abs(xq(t))) - 1;
%     comp = 2*(real(x(t)*exp(1j*pi/4)) < imag(x(t)*exp(1j*pi/4))) - 1;
    
    y(t) = xid(t)*xqd(t)*comp;

%     y(t) = xqd(t)*xi(t) - xid(t)*xq(t);

%     y(t) = xi.*xq;

    yf(t) = sum(numz.*y(t:-1:t-numLen+1)) - sum(yf(t-1:-1:t-denLen+1).*denz(2:end));
end

figure
subplot(211)
plot(y)
subplot(212)
plot(yf)


scatterplot(xn)
title('input')
scatterplot(xr)
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
    
    
    


