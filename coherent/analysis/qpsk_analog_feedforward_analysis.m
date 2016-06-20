%% Analysis of carrier phase recovery for QPSK using feedforward technique

clear, clc, close all

addpath ../../f/

Nsymb = 2^10;
Mct = 15;
N = Nsymb*Mct;
Rs = 56e9;
fs = Rs*Mct;
ts = 0:1/fs:(N-1)*1/fs;
f = -fs/2:fs/N:fs/2-fs/N;

Laser = laser(1310e-9, 0, -150, 20000e3);

phase_noise = true;
awgn_noise = ~true;
frequency_offset = ~true;

data = randi([1 4], [1 Nsymb]);

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
    phin = 0;
    xn = x;
end
xn = x.*exp(1j*phin);

varN = 0.1;
if awgn_noise
    xn = xn + sqrt(varN/2)*randn(size(xn)) + 1j*sqrt(varN/2)*randn(size(xn));
end

foff = 1e9;
if frequency_offset   
    xn = freqshift(xn, 0:1/fs:(length(xn)-1)/fs, foff);
end

xr = zeros(size(x));
y = zeros(size(x));
yf = zeros(size(x));

nand = @(a,b) -(2*(a > 0 && b > 0) - 1);

xn4d = sign(imag(xn.^4));
q = zeros(size(x))+1;
qbar = zeros(size(x))+1;
for t = Mct+1:length(x)     
    q(t) = nand(nand(xn4d(t),qbar(t-4)), qbar(t-2));
    qbar(t) = nand(nand(xn4d(t), q(t-4)), q(t-2));
end
    
 
    
    
    
   
%     xi(t) = real(xn(t));
%     xq(t) = imag(xn(t));
%     
%     xid(t) = sign(xi(t));
%     xqd(t) = sign(xq(t));
%     comp = 2*(abs(xi(t)) < abs(xq(t))) - 1;
% %     comp = 2*(real(x(t)*exp(1j*pi/4)) < imag(x(t)*exp(1j*pi/4))) - 1;
%     
%     y(t) = xid(t)*xqd(t)*comp;
% 
% %     y(t) = xqd(t)*xi(t) - xid(t)*xq(t);
% 
% %     y(t) = xi.*xq;


%     comp = 2*(real(x(t)*exp(1j*pi/4)) < imag(x(t)*exp(1j*pi/4))) - 1;
%     
%     tmp = not(xor(xid(t), xqd(t)));
%     y(t) = not(xor(tmp, comp));
%     y(t) = 2*y(t) - 1;
   
%     y(t) = xqd(t)*xi(t) - xid(t)*xq(t);    


figure
hold on
plot(xn4d, '-r')
% plot(q)
% plot(sqrt(2)*sin(phin), '--k')
% subplot(212)
% plot(yf)


% scatterplot(xn)
% title('input')
% scatterplot(xr)
% title('output')

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
    
    
    


