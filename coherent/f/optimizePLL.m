function [wnOpt, minVarPhiError, SNRdB] = optimizePLL(csi, Delay, noiseParam, Ncpr, sim, verbose)
%% Optimizes PLL natural frequency given csi (damping), Delay, linewidth for a target BER
% Inputs:
% - csi: loop filter damping factor
% - Delay: total loop delay
% - noiseParam: [total linewidth i.e., transmitter laser + LO,
% flickerNoise/f^3 is the one-sided phase noise PSD due to flicker noise (optional)]
% - Ncpr: number of polarizations used in carrier recovery
% - sim: simulation struct
% - verbose: plot
% Outputs:
% - wnOpt: optimal natural frequency
% - minVarPhiError: value of phase error variance at optimal wn
% - SNRdB: SNR in dB at each optimization was computed

wnmax = min(csi/Delay, 2*pi*1e9); % rule of thumb for the relaxation frequency that leads to maximum phase error variance
M = sim.ModFormat.M;
BER = sim.BERtarget;
SNRdB = SNRreq(BER, M, sim.ModFormat.type);

[wnMradpsOpt, minVarPhiError, exitflag] = fminbnd(@(wnMradps) ...
    phase_error_variance(csi, wnMradps*1e6, Ncpr, Delay, noiseParam, SNRdB, sim.ModFormat.Rs),...
    0, wnmax/1e6);

wnOpt = wnMradpsOpt*1e6;

if exitflag ~= 1
    warning('PLL relaxation frequency optimization ended with exitflag = %d\n', exitflag);
end

if exist('verbose', 'var') && verbose
    wn = linspace(0, wnmax);
    varPhiError = phase_error_variance(csi, wn, Ncpr, Delay, noiseParam, SNRdB, sim.ModFormat.Rs);
    [minVarPhiError, nPN, nAWGN, nFlicker] = phase_error_variance(csi, wnOpt, Ncpr, Delay, noiseParam, SNRdB, sim.ModFormat.Rs);
    
    fprintf('Contribution on phase error variance at SNRdB = %.2f:\nPN/AWGN = %.3f\nFlicker/AWGN = %.3f\n', SNRdB,...
        nPN/nAWGN, nFlicker/nAWGN);
    
    figure(300), hold on
    h = plot(wn/1e9, varPhiError);
    plot(wnOpt/1e9, minVarPhiError, 'o', 'Color', get(h, 'Color'))
    xlabel('\omega_n (Grad/s)')
    ylabel('Phase error variance (rad^2)')
end