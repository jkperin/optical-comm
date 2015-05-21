%% Calculate optimal Gain for APD

% Input:
% mpam : signal struct (bw = bandwidth, M = number of levels)
% tx : transmitter struct (Prec = average received power, RIN = RIN in dB/Hz)
% rx : receiver struct (Gapd = apd gain, ka = avalanche efficiency, R = 
% responsivity)
% sim : simulation paramters (struct) (SERtarget = target SER)
% maxGain: maximum Gain (may be limited by the Gain x Bandwidth product or
% not)
% levelspacing: {'uniform', 'nonuniform'}

% Output:
% Gapd_opt = optimal APD gain

function Gapd_opt = calcOptAPDGain(mpam, tx, apd, sim, maxGain)

if strcmp(sim.levelSpacing, 'uniform')
    % Optmize gain for uniform spacing: find Gapd that minimizes the required
    % average power (Prec) to achieve a certain target SER.
    [Gapd_opt, ~, exitflag] = fminbnd(@(Gapd) fzero(@(PrecdBm) ber_mpam(mpam, tx, apd, sim, Gapd, 1e-3*10^(PrecdBm/10)) - sim.BERtarget, -20), 1, maxGain);

elseif strcmp(sim.levelSpacing, 'nonuniform')
    % Optimal level spacing
    [Gapd_opt, ~, exitflag] = fminbnd(@(Gapd) mean(calcOptLevelSpacing(mpam, tx, apd, sim, Gapd)), 1, maxGain);

else
    error('Invalid Option')
end
    

if exitflag ~= 1
    warning(sprintf('APD gain optimization did not converge (exitflag = %d)\n', exitflag))
end    