clear, clc, close all  

M = 4;
SNRsdB = 10; % dB
Ts = 1/56e9;
LineWidth = 1e6;

SNRs = 10^(SNRsdB/10);

const = qammod(0:M-1, M, 0, 'Gray');
eta = 1/2*mean(abs(const).^2)*mean(abs(1./const).^2);

varPN = 2*pi*LineWidth*Ts;
varNp = eta*1/SNRs;

r = varPN/varNp;
alpha = (1 + r/2) - sqrt((1 + r/2)^2-1);

n = -20:20;
wn = alpha*r/(1 - alpha^2)*alpha.^(abs(n));

stem(n, wn/max(wn))
