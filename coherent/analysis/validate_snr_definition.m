%% Validate BER calculation
clear, clc, close all

N = 2^14;

Mpam = 2;
Mqam = 4;
dataTXpam = randi([0 Mpam-1], [N 1]);
dataTXqam = randi([0 Mqam-1], [N 1]);

EbN0dB = 0:10;
EbN0 = 10.^(EbN0dB/10);

for k = 1:length(EbN0dB)
    xpam = pammod(dataTXpam, Mpam, 0, 'gray'); % EsRs = 1, EbRs = 1
    xqam = qammod(dataTXqam, Mqam, 0, 'gray'); % EsRs = 2, EbRs = 1
    
    npam = sqrt(0.5/EbN0(k))*randn(N, 1);
    nqam = sqrt(0.5/EbN0(k))*randn(N, 1) + 1j*sqrt(0.5/EbN0(k))*randn(N, 1);
        
    ypam = xpam + npam;
    yqam = xqam + nqam;
    
    dataRXpam = pamdemod(ypam, Mpam, 0, 'gray');
    dataRXqam = qamdemod(yqam, Mqam, 0, 'gray');
    
    [~, BERpam(k)] = biterr(dataTXpam, dataRXpam);
    [~, BERqam(k)] = biterr(dataTXqam, dataRXqam);
end

figure, hold on
plot(EbN0dB, log10(BERpam), '-ob', EbN0dB, log10(BERqam), '-or')
plot(EbN0dB, log10(berawgn(EbN0dB, 'pam', Mpam)), '-b',...
    EbN0dB, log10(berawgn(EbN0dB, 'qam', Mqam)), '-r')


figure, hold on
SNRdBpam = EbN0dB + 10*log10(log2(Mpam)) + 10*log10(2);
SNRdBqam = EbN0dB + 10*log10(log2(Mqam));
plot(SNRdBpam, log10(BERpam), '-ob', SNRdBqam, log10(BERqam), '-or')
plot(SNRdBpam, log10(berawgn(EbN0dB, 'pam', Mpam)), '-b',...
    SNRdBqam, log10(berawgn(EbN0dB, 'qam', Mqam)), '-r')