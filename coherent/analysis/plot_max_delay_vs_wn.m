%% Maximum delay vs loop filter natural frequency
clear, clc, close all

addpath ../f/
addpath ../../f/

wn = 2*pi*1e6*linspace(20, 180, 21);
totalLinewidth = 1e3*[50, 100, 200, 400];

sim.ModFormat = QAM(4, 2*112e9);
sim.BERtarget = 1.8e-4;
SNRdBpen = 0.5;
csi = sqrt(2)/2;
Ncpr = 2;
SNRdBref = SNRreq(sim.BERtarget, sim.ModFormat.M, 'QAM'); % Required SNR to achieve target BER

for k = 1:length(totalLinewidth)
    parfor i = 1:length(wn)
        fprintf('l = %.2f kHz, wn = %.2f MHz\n', totalLinewidth(k)/1e3, wn(i)/1e6);
        maxDelay(k, i) = max_loop_delay(SNRdBref + SNRdBpen, csi, totalLinewidth(k), Ncpr, sim, wn(i), false);
        
        [~, varPN, varAWGN] = phase_error_variance(csi, wn(i), Ncpr, maxDelay(k, i), totalLinewidth(k), SNRdBref + SNRdBpen, sim.ModFormat.Rs, false);
        fprintf('> varPN/varAWGN = %.2f\n', varPN/varAWGN)
    end
    % last one is at the optimal wn
    [maxDelay(k, length(wn) + 1), wnopt(k)] = max_loop_delay(SNRdBref +SNRdBpen, csi, totalLinewidth(k), Ncpr, sim);
end

figure, hold on, box on
plot(wn/(2*pi*1e6), 1e12*maxDelay(:, 1:end-1))
plot(wnopt/(2*pi*1e6), 1e12*maxDelay(:, end), 'o')
xlabel('Natural frequency (MHz)')
ylabel('Maximum delay for 0.5-dB SNR penalty (ps)')
axis([wn([1 end])/(2*pi*1e6) 0 4000])
legend('linewidth = 50 kHz', '100', '200', '400')