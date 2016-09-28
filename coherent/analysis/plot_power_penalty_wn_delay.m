%% Plot receiver sensitivity penalty due to non-optimal loop bandwidth and delay
% clear, clc, close all
clear, clc, close all

addpath ../../f/
addpath ../f/

Qpsk = QAM(4, 2*112e9);
BERtarget = 1.8e-4;
sim.BERtarget = BERtarget;
sim.ModFormat = Qpsk;
SNRdBref = SNRreq(BERtarget, Qpsk.M, 'QAM');
R = 1;
q = 1.60217662e-19;

csi = sqrt(2)/2;
Ncpr = 2;

PrxdBmref = 10*log10(10^(SNRdBref/10)*2*q*Qpsk.Rs/(R*1e-3))

totalLinewidth = 2e3*(100:100:1e3);

Delays = (0:100:1000)*1e-12;

for l = 1:length(totalLinewidth)
    for kk = 1:length(Delays)
        [wnOpt, minVarPhiError, SNRdB] = optimizePLL(csi, Delays(kk), totalLinewidth(l), Ncpr, sim);
        
        [SNRdBpen(l, kk), ~, exitflag] = fzero(@(SNRdB)...
            log10(ber_qpsk_imperfect_cpr(SNRdB, phase_error_variance(csi, wnOpt, Ncpr, Delays(kk), totalLinewidth(l), SNRdB, Qpsk.Rs))) - log10(BERtarget), SNRdBref);

        if exitflag ~= 1
            warning('SNRdB calculation finished with exitflag = %d', exitflag);
        end

        PrxreqdBm(l, kk) = 10*log10(10^(SNRdBpen(l, kk)/10)*2*q*Qpsk.Rs/1e-3);
    end
end

PpendB = PrxreqdBm - PrxdBmref;
Pmax = 0.5;

figure, hold on
for l = 1:length(totalLinewidth)
    p = polyfit(1e12*Delays, PpendB(l, :), 2);
    plot(1e12*Delays, PpendB(l, :));
    plot(1e12*Delays, polyval(p, 1e12*Delays), '--');
    drawnow
    
    p(end) = p(end) - Pmax;
    r = roots(p)
    
    maxDelay(l) = min(r(r > 0));
end

figure, plot(totalLinewidth, maxDelay)

