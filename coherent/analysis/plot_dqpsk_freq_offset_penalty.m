%% Plot DPSK frequency offset penalty

clear, clc, close all

addpath ../f/
addpath ../../f

Nsymb = 2^16;
Rb = 112e9;
ros = 1;
N = Nsymb*ros;
BERtarget = 1.8e-4;

pulse_shape = select_pulse_shape('rect', ros);
Dpsk = DPSK(4, Rb, pulse_shape);

fs = Dpsk.Rs*ros;
[f, t] = freq_time(N, fs);

foff = linspace(0, 5, 100)*1e9;

SNRdBreq = SNRreq(BERtarget, Dpsk.M, Dpsk.type) % required SNR for target BER

SNRdBnew =zeros(size(foff));
for k = 1:length(foff)
    [SNRdBnew(k), ~, exitflag] = fzero(@(SNRdB) log10(Dpsk.ber_freq_offset(SNRdB, foff(k))) - log10(BERtarget), SNRdBreq);
    
    if exitflag ~= 1
        exitflag 
        SNRdBnew(k) = 0;
    end
end

SNRdBpen = SNRdBnew - SNRdBreq;
R = 1; % photodiodes responsivity
q = 1.60217662e-19; % electron charge
SNRdB2PrxdBm = @(SNRdB) 10*log10(10.^(SNRdB/10)*2*q*Dpsk.Rs/(R*1e-3));
Ppen = SNRdB2PrxdBm(SNRdBnew) - SNRdB2PrxdBm(SNRdBreq);

figure(1), hold on, box on
plot(foff/1e9, Ppen, '-k', 'LineWidth', 2)
xlabel('Frequency offset (GHz)', 'FontSize', 12)
ylabel('Power penalty (dB)', 'FontSize', 12)
grid on
axis([foff([1 end])/1e9 0 10])
m2tikz = matlab2tikz(gca);
m2tikz.write('dqpsk_freq_offset_penalty.tikz');