%% Validate mutual_info.m
clear, clc, close all

addpath ../mpam/

N = 2^15;
SNRdB = 0:20;
SNR = 10.^(SNRdB/10);

MPAM = PAM(4, 1);
MPAM = MPAM.unbias;
Xpam = [MPAM.mod(randi([0 MPAM.M-1], [N, 1])).',...
    MPAM.mod(randi([0 MPAM.M-1], [N, 1])).'];
Pspam = MPAM.Ppam;

MQAM = QAM(4, 1);
Xqam = MQAM.mod(randi([0 MQAM.M-1], [N, 1]));
Xqam = [real(Xqam) imag(Xqam)];
Psqam = MQAM.Pqam;

for k = 1:length(SNR)
    % Parallel M-PAM
    Ypam = Xpam + sqrt(Pspam/(SNR(k)))*randn(N, 2);
    MIpam(k) = mutual_info(Xpam, Ypam, 10);
    
    % M-QAM
    Yqam = Xqam + sqrt(Psqam/(2*SNR(k)))*randn(N, 2);
    MIqam(k) = mutual_info(Xqam, Yqam, 10);
end

figure, box on, hold on
plot(SNRdB, MIpam, 'LineWidth', 2, 'DisplayName', sprintf('2 x %d-PAM', MPAM.M))
plot(SNRdB, MIqam, 'LineWidth', 2, 'DisplayName', sprintf('%d-QAM', MQAM.M))
plot(SNRdB, 0.5*log2(1 + SNR), 'LineWidth', 2, 'DisplayName', 'Shannon 1 DOF')
plot(SNRdB, log2(1 + SNR), 'LineWidth', 2, 'DisplayName', 'Shannon 2 DOF')
xlabel('SNR (dB)')
ylabel('Spectral efficiency (bits/s/Hz)')
legend('-dynamiclegend')
