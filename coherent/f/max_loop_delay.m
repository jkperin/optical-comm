function [maxDelay, wnSet] = max_loop_delay(SNRdB, csi, noiseParam, Ncpr, sim, wnOffFactor, optimizeWn)
%% Calculate maximum loop delay to achieve a given target BER with SNR = SNRdBpen
% Inputs:
% - SNRdB: SNR in dB including maximum accepted penalty due to phase
% error
% - csi: Loop filter damping constant
% - noiseParam: [total linewidth i.e., transmitter laser + LO,
% flickerNoise/f^3 is the one-sided phase noise PSD due to flicker noise (optional)]
% - Ncpr: number of polarizations used in CPR {1, 2}
% - sim: simulation struct containing modulation format and target BER
% - wnOffFactor (optional): in calculating penalty wn = wnOffFactor*wnOpt
% Calculate new SNR necessary to make system operate at target BER
% taking into account errors due to imperfect carrier phase
% recovery
% - optimizeWn (optinal, default = true): whether to use optimal wn or the value provided in wnOffFactor

if not(exist('wnOffFactor', 'var'))
    wnOffFactor = 1; % i.e., wn = optimal wn
end

if not(exist('optimizeWn', 'var'))
    optimizeWn = true;
end

[maxDelay, ~, exitflag] = fzero(@(delay) log10(calc_ber(abs(delay)*1e-12)) - log10(sim.BERtarget), 1000);

if exitflag ~= 1
    warning('max_loop_delay: optimization exited with exitflag = %d', exitflag);
    maxDelay = NaN;
    wnOpt = NaN;
else
    maxDelay = 1e-12*abs(maxDelay);
    wnOpt = optimizePLL(csi, maxDelay, noiseParam, Ncpr, sim);
end

    function ber = calc_ber(delay)
        if optimizeWn
            wnOpt = optimizePLL(csi, delay, noiseParam, Ncpr, sim);
            wnSet = wnOffFactor*wnOpt;
        else
            wnSet = wnOffFactor; % in this case wnOffFactor is actually the chosen wn
        end
            
        ber = ber_qpsk_imperfect_cpr(SNRdB,...
            phase_error_variance(csi, wnSet, Ncpr, delay, noiseParam, SNRdB, sim.ModFormat.Rs, false));
    end
end

