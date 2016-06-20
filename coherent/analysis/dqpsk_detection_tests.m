%% Detection of DQPSK using 1-bit ADC

clear, clc, close all

addpath ../../f/

%
Nsymb = 2^14;
Mct = 1;
N = Nsymb*Mct;
Rs = 56e9;
fs = Rs*Mct;

Laser = laser(1310e-9, 0, -150, 200e3);

phase_noise = true;
awgn_noise = true;
frequency_offset = ~true;

constellation = pi/4*[1 3 5 7];
data = randi([1 4], [Nsymb 1]);

% xi = cos(constellation(data));
% xq = sin(constellation(data));
pchange = [0 -pi/2 pi/2 pi];
phaseChangeIdeal = pchange(data);
x = zeros(1, Nsymb);
x(1) = exp(1j*pi/4);
for t = 2:length(x)
    x(t) = x(t-1)*exp(1j*phaseChangeIdeal(t));
end

xi = real(x);
xq = imag(x);

% xi = reshape(repmat(xi, Mct, 1), [], N);
% xq = reshape(repmat(xq, Mct, 1), [], N);

x = xi + 1j*xq;

if phase_noise
    [xn, phin] = Laser.addPhaseNosie(x, fs);
else
    xn = x;
end

varN = 0.001;
if awgn_noise
    xn = xn + sqrt(varN/2)*randn(size(xn)) + 1j*sqrt(varN/2)*randn(size(xn));
end

foff = 1e9;
if frequency_offset   
    xn = freqshift(xn, 0:1/fs:(length(xn)-1)/fs, foff);
end

%% Detection
xr = xn;
xid = real(xr);
xqd = imag(xr);
for t = Mct+1:length(x)   
    xi(t) = real(xr(t));
    xq(t) = imag(xr(t));
    
    xid(t) = xi(t) >= 0;
    xqd(t) = xq(t) >= 0;

    si = xor(xid(t), xid(t-Mct));
    sq = xor(xqd(t), xqd(t-Mct));
    a = xid(t-Mct);
    b = xqd(t-Mct);
    c = xid(t);
    d = xqd(t);
    ctrl = ((~ a) && b && (~ c) && (~ d))...
        || (a && (~ b) && (~ c) && (~ d))...
        || ((~ a) && b && c && d)...
        || (a && (~ b) && c && d);
    di(t) = xor(si, ctrl);
    dq(t) = xor(sq, ctrl);
    1;
end

dataRX = zeros(size(x))+1;
dataRX(di  == 0 & dq == 1) = 2;
dataRX(di  == 1 & dq == 0) = 3;
dataRX(di  == 1 & dq == 1) = 4;


[~, ber] = biterr(data(512:end).', dataRX(512:end))


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
    
    
    


