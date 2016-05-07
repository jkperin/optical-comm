clear, clc, close all

addpath ../f/

Nsymb = 2^10;
Mct = 15;
N = Nsymb*Mct;
Rs = 56e9;
fs = Rs*Mct;

Laser = laser(1310e-9, 0, -150, 400e3);

phase_noise = true;
awgn_noise = true;
frequency_offset = true;

% constellation = pi/4*[1 3 5 7];
constellation = pi*[0 1];
data = randi([1 2], [Nsymb 1]);

xi = cos(constellation(data));
xq = sin(constellation(data));

xi = reshape(repmat(xi, Mct, 1), [], N);
xq = reshape(repmat(xq, Mct, 1), [], N);

x = xi + 1j*xq;

if phase_noise
    [xn, phin] = Laser.addPhaseNosie(x, fs);
else
    xn = x;
end

if awgn_noise
    xn = xn + sqrt(0.01/2)*randn(size(xn)) + 1j*sqrt(0.01/2)*randn(size(xn));
end

if frequency_offset    
    xn = freqshift(xn, 0:1/fs:(length(xn)-1)/fs, 10e9);
end

% Loop filter
csi = 1/sqrt(2);                                                    % damping coefficient of second-order loop filter
wn = 2*pi*2.4e9;                                                    % relaxation frequency of second-order loop filter: optimized using optimize_PLL.m
nums = [2*csi*wn wn^2];
dens = [1 0 0]; % descending powers of s

[numz, denz] = impinvar(nums, dens, fs);
numLen = length(numz);
denLen = length(denz);

xr = zeros(size(x));
y = zeros(size(x));
yf = zeros(size(x));
for t = numLen+1:length(x)
    xr(t) = exp(1j*yf(t-1))*xn(t);
    
    xi = real(xr(t));
    xq = imag(xr(t));
    
    xid(t) = (xi > 0);
    xqd(t) = (xq > 0);

    y(t) = xor(xid(t), xqd(t));
    y(t) = 2*y(t) - 1;

    yf(t) = sum(numz.*y(t:-1:t-numLen+1)) - sum(yf(t-1:-1:t-denLen+1).*denz(2:end));
end

figure
plot(yf)

scatterplot(xn)
title('input')
scatterplot(xr)
title('output')


