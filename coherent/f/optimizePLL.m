function [wnOpt, wn, phiError, SNRdB] = optimizePLL(csi, Delay, totalLinewidth, Ncpr, sim, verbose)
%% Optimizes PLL relaxation frequency given csi (damping), Delay, linewidth for a target BER

Ts = 1/sim.ModFormat.Rs;
M = sim.ModFormat.M;
BER = sim.BERtarget;
SNRdB = SNRreq(BER, M, sim.ModFormat.type);
SNR = 10^(SNRdB/10);

wn = 2*pi*1e9*(0:0.01:0.3); % at most ~ 1.8 Grad/s
w = 2*pi*linspace(-0.5/Ts, 0.5/Ts, 2^12);
phiError = zeros(size(wn));
for k = 1:length(wn)
    numFs = [2*csi*wn(k) wn(k)^2];
    denFs = [1 0];
    
    Fw = polyval(numFs, 1j*w)./polyval(denFs, 1j*w);
    H1 = 1j*w + exp(-1j*w*Delay).*Fw;
    
    phiError(k) = totalLinewidth*trapz(w, 1./abs(H1).^2) + Ts/(2*pi*Ncpr*2*SNR)*trapz(w, abs(Fw./H1).^2);
end

[minPhiError, ind] = min(phiError);
wnOpt = wn(ind);

if exist('verbose', 'var') && verbose
    numFs = [2*csi*wnOpt wnOpt^2];
    denFs = [1 0];
    Fw = polyval(numFs, 1j*w)./polyval(denFs, 1j*w);
    H1 = 1j*w + exp(-1j*w*Delay).*Fw;
    PNoverAWGN = totalLinewidth*trapz(w, 1./abs(H1).^2)/(Ts/(2*pi*Ncpr*2*SNR)*trapz(w, abs(Fw./H1).^2));    
    fprintf('Contribution of phase noise vs AWGN on PLL phase error at receiver sensitivity (SNRdB = %.2f):\nPN/AWGN = %.3f\n', SNRdB, PNoverAWGN);
    
    figure(300), hold on
    h = plot(wn/1e9, phiError);
    plot(wnOpt/1e9, minPhiError, 'o', 'Color', get(h, 'Color'))
    xlabel('\omega_n (Grad/s)')
    ylabel('Phase error variance (rad^2)')
end