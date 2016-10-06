function [maxDelay, wnOpt] = max_loop_delay(SNRdBpen, csi, totalLinewidth, Ncpr, sim, wnOffFactor)
%% Calculate maximum loop delay to achieve a given target BER with SNR = SNRdBpen
% Inputs:
% - SNRdBpen: SNR in dB including maximum accepted penalty due to phase
% error
% - csi: Loop filter damping constant
% - totalLinewidth: total linewidth i.e., combined of transmitter laser and
% LO
% - Ncpr: number of polarizations used in CPR {1, 2}
% - sim: simulation struct containing modulation format and target BER
% - wnOffFactor (optional): in calculating penalty wn = wnOffFactor*wnOpt
% Calculate new SNR necessary to make system operate at target BER
% taking into account errors due to imperfect carrier phase
% recovery

if not(exist('wnOffFactor', 'var'))
    wnOffFactor = 1; % i.e., wn = optimal wn
end

[maxDelay, ~, exitflag] = fzero(@(delay) log10(calc_ber(abs(delay)*1e-12)) - log10(sim.BERtarget), 100);

if exitflag ~= 1
    warning('max_loop_delay: optimization exited with exitflag = %d', exitflag);
    maxDelay = NaN;
    wnOpt = NaN;
else
    maxDelay = 1e-12*abs(maxDelay);
    wnOpt = optimizePLL(csi, maxDelay, totalLinewidth, Ncpr, sim);
end

    function ber = calc_ber(delay)
        wnOpt = optimizePLL(csi, delay, totalLinewidth, Ncpr, sim);
        ber = ber_qpsk_imperfect_cpr(SNRdBpen,...
            phase_error_variance(csi, wnOffFactor*wnOpt, Ncpr, delay, totalLinewidth, SNRdBpen, sim.ModFormat.Rs, false));
    end
end

