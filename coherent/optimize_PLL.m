clear, clc, close all
addpath f/

Dnu = 0.2e6;
Ts = 1/56e9;
csi = 1/sqrt(2);
M = 4;
BER = 1e-4;
SNRdB = SNRreq(BER, M, 'QAM')
SNR = 10^(SNRdB/10);
etac = mean(abs(qammod(0:M-1, M)).^2)/mean(1./abs(qammod(0:M-1, M)).^2);

varPN = 2*pi*Dnu*Ts;

wn = 1e9*(0:0.05:1);
w = linspace(-0.5/Ts, 0.5/Ts, 1e3);
s = tf('s');
for k = 1:length(wn)
    Fs = 2*csi*wn(k) + wn(k)^2/s;
    
    [num1, den1] = tfdata(1/(s + Fs));
    num1 = cell2mat(num1);
    den1 = cell2mat(den1);
    H1 = freqs(num1, den1, w);
%     figure, plot(w/1e9, abs(H1).^2)
    GammaPN(k) = 2*csi*wn(k)/pi*trapz(w, abs(H1).^2);
    
    [num2, den2] = tfdata(Fs/(s + Fs));
    num2 = cell2mat(num2);
    den2 = cell2mat(den2);
    H2 = freqs(num2, den2, w);
%     figure, plot(w/1e9, abs(H2).^2)
    GammaAWGN(k) = 2*csi/(pi*(1 + 4*csi^2)*wn(k))*trapz(w, abs(H2).^2);
    
    phiError(k) = varPN/(2*csi*wn(k)*Ts)*GammaPN(k) + (1 + 4*csi^2)*wn(k)*Ts/(4*csi)*(etac/(2*SNR))*GammaAWGN(k);
end

figure, hold on
plot(wn/1e9, phiError)

[minPhiError, ind] = min(phiError);
wnOpt = wn(ind)
plot(wnOpt/1e9, minPhiError, 'o')
xlabel('\omega_n (GHz)')
ylabel('Phase error variance')

figure
G = Fs/s;
G = G/(1 + G);
step(G)
