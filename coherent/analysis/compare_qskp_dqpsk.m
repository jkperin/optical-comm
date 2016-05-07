    %% Compare QPSK and DPSK in AWGN
clear, clc, close all

M = 4;
N = 2^12;

x = randi([0 M-1], [1 N]);
yQPSK = qammod(x, M, 0, 'Gray');
yDPSK = sqrt(2)*exp(1j*pi/4)*dpskmod(x, M, 0, 'Gray');

SNRdB = 6:20;
SNR = 10.^(SNRdB/10);
validInd = 1:N;
EsQAM = mean(abs(qammod(0:M-1, M)).^2);
EsDPSK = mean(abs(sqrt(2)*exp(1j*pi/4)*dpskmod(x, M)).^2);
% SNR = Es/N0, N0 = 2sigma^2

for k = 1:length(SNRdB)
    yQPSKn = yQPSK + sqrt(EsQAM/(SNR(k)))*(randn(1, N) + 1j*randn(1, N));
    yDPSKn = yDPSK + sqrt(EsDPSK/(SNR(k)))*(randn(1, N) + 1j*randn(1, N));
    
    xhatQPSK = qamdemod(yQPSKn, M, 0, 'gray');
    xhatDPSK = dpskdemod(1/sqrt(2)*exp(-1j*pi/4)*yDPSKn(validInd), M, 0, 'gray').';
    
    [~, berQPSK(k)] = biterr(xhatQPSK, x);
    [~, berDPSK(k)] = biterr(xhatDPSK, x(validInd));
end

figure, hold on, box on
plot(SNRdB, log10(berQPSK), '-o')
plot(SNRdB, log10(berDPSK), '-o')
plot(SNRdB, log10(berawgn(SNRdB-3-10*log10(log2(M)), 'qam', M)), 'k')
plot(SNRdB, log10(berawgn(SNRdB-3-10*log10(log2(M)), 'dpsk', M)), '--k')
axis([SNRdB([1 end]) -8 0])