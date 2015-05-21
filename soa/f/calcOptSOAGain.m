function Gsoa_opt = calcOptSOAGain(mpam, tx, soa, rx, sim, maxGain)

if strcmp(sim.levelSpacing, 'uniform')
    % Optmize gain for uniform spacing: find Gapd that minimizes the required
    % average power (Prec) to achieve a certain target SER.
    [Gsoa_opt, ~, exitflag] = fminbnd(@(Gsoa) fzero(@(PrecdBm) mpam_ber_soa(mpam, soa, rx, Gsoa, PrecdBm) - sim.BERtarget, -20), 1, maxGain);

elseif strcmp(sim.levelSpacing, 'nonuniform')
    % Optimal level spacing
    [Gsoa_opt, ~, exitflag] = fminbnd(@(Gsoa) mean(calcOptLevelSpacingAmpl(mpam, tx, soa, rx, sim, Gsoa)), 1, maxGain);

else
    error('Invalid Option')
end
    

if exitflag ~= 1
    warning(sprintf('APD gain optimization did not converge (exitflag = %d)\n', exitflag))
end    