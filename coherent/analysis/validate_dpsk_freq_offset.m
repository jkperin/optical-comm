%% Impact of frequency offset in DPSK
clear, clc, close all

addpath ../../f

Nsymb = 2^16;
Rb = 112e9;
ros = 1;
N = Nsymb*ros;

pulse_shape = select_pulse_shape('rect', ros);
Dpsk = DPSK(4, Rb, pulse_shape);

fs = Dpsk.Rs*ros;
[f, t] = freq_time(N, fs);

foff = 2e9;

SNRdB = 5:20;
SNR = 10.^(SNRdB/10);

dataTX = randi([0 Dpsk.M-1], [1 Nsymb]);

for k = 1:length(SNRdB)
    x = Dpsk.signal(dataTX);
    Ps = mean(abs(x).^2);
    Pn = Ps/SNR(k);
    n = sqrt(Pn/2)*(randn(size(x)) + 1j*randn(size(x)));
    
    y = x.*exp(1j*2*pi*foff*t) + n;
    
    dataRX = Dpsk.demod(y);
    
    [~, ber(k)] = biterr(dataTX, dataRX);
end

bertheory = Dpsk.ber_freq_offset(SNRdB, foff);

figure(1), hold on, box on
plot(SNRdB, log10(ber), '-o')
plot(SNRdB, log10(berawgn(SNRdB - 10*log10(log2(Dpsk.M)), 'dpsk', Dpsk.M)));
plot(SNRdB, log10(bertheory), '--k');
xlabel('SNR (dB)')
ylabel('log_{10} (BER)')
legend('Counted', 'AWGN', 'Theory freq offset')
axis([SNRdB([1 end]) -8 0])

    
