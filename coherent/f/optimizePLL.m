function [wnOpt, minVarPhiError, SNRdB] = optimizePLL(csi, Delay, totalLinewidth, Ncpr, sim, verbose)
%% Optimizes PLL relaxation frequency given csi (damping), Delay, linewidth for a target BER

wnmax = min(csi/Delay, 2*pi*1e9); % rule of thumb for the relaxation frequency that leads to maximum phase error variance
M = sim.ModFormat.M;
BER = sim.BERtarget;
SNRdB = SNRreq(BER, M, sim.ModFormat.type);

[wnMradpsOpt, minVarPhiError, exitflag] = fminbnd(@(wnMradps) ...
    phase_error_variance(csi, wnMradps*1e6, Ncpr, Delay, totalLinewidth, SNRdB, sim.ModFormat.Rs),...
    0, wnmax/1e6);

wnOpt = wnMradpsOpt*1e6;

if exitflag ~= 1
    warning('PLL relaxation frequency optimization ended with exitflag = %d\n', exitflag);
end

if exist('verbose', 'var') && verbose
    wn = linspace(0, wnmax);
    varPhiError = phase_error_variance(csi, wn, Ncpr, Delay, totalLinewidth, SNRdB, sim.ModFormat.Rs);
    [minVarPhiError, nPN, nAWGN] = phase_error_variance(csi, wnOpt, Ncpr, Delay, totalLinewidth, SNRdB, sim.ModFormat.Rs);
    PNoverAWGN = nPN/nAWGN;
    
    fprintf('Contribution of phase noise vs AWGN on PLL phase error at receiver sensitivity (SNRdB = %.2f):\nPN/AWGN = %.3f\n', SNRdB, PNoverAWGN);
    
    figure(300), hold on
    h = plot(wn/1e9, varPhiError);
    plot(wnOpt/1e9, minVarPhiError, 'o', 'Color', get(h, 'Color'))
    xlabel('\omega_n (Grad/s)')
    ylabel('Phase error variance (rad^2)')
end