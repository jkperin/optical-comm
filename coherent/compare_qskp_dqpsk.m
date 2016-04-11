    %% Compare QPSK and DQPSK in AWGN
clear, clc, close all

M = 4;
N = 2^12;

x = randi([0 M-1], [1 N]);
yQPSK = qammod(x, M, 0, 'Gray');
yDQPSK = sqrt(2)*exp(1j*pi/4)*dpskmod(x, M, 0, 'Gray');

SNRdB = 0:15;
SNR = 10.^(SNRdB/10);
validInd = 324:N;

for k = 1:length(SNRdB)
    yQPSKn = yQPSK + sqrt(0.5/SNR(k))*(randn(1, N) + 1j*randn(1, N));
    yDQPSKn = yDQPSK + sqrt(0.5/SNR(k))*(randn(1, N) + 1j*randn(1, N));
    
    xhatQPSK = qamdemod(yQPSKn, M, 0, 'gray');
    xhatDQPSK = dpskdemod(1/sqrt(2)*exp(-1j*pi/4)*yDQPSKn(validInd), M, 0, 'gray').';
    
    [~, berQPSK(k)] = biterr(xhatQPSK, x);
    [~, berDQPSK(k)] = biterr(xhatDQPSK, x(validInd));
end

figure, hold on, box on
plot(SNRdB, log10(berQPSK), '-o')
plot(SNRdB, log10(berDQPSK), '-o')
plot(SNRdB, log10(berawgn(SNRdB-log10(log2(M)), 'qam', M)), 'k')
plot(SNRdB, log10(berawgn(SNRdB-log10(log2(M)), 'dpsk', M)), '--k')
axis([SNRdB([1 end]) -8 0])