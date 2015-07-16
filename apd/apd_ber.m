%% Calculate BER of unamplified IM-DD system with APD detector 
% BER is calculated via montecarlo simulation and through saddlepoint approximation
% to calculate tail probabilities given the moment generating funcion.

function [ber, mpam] = apd_ber(mpam, tx, fiber, apd, rx, sim)

dBm2Watt = @(x) 1e-3*10.^(x/10);

% Level Spacing Optimization
if strcmp(mpam.level_spacing, 'optimized')
    % Noises variance
    varTherm = rx.N0*rx.elefilt.noisebw(sim.fs)/2; % variance of thermal noise

    % Variance of RIN
    if sim.RIN
        var_rin = @(P) 10^(tx.RIN/10)*P.^2*rx.elefilt.noisebw(sim.fs);
    else
        var_rin = @(P) 0;
    end

    % Shot noise variance = Agrawal 4.4.17 (4th edition)
    calc_noise_std = @(Plevel) sqrt(varTherm + var_rin(Plevel/apd.Gain) + apd.var_shot(Plevel/apd.Gain, rx.elefilt.noisebw(sim.fs)/2));

    % Optimize levels using Gaussian approximation
    mpam.optimize_level_spacing_gauss_approx(sim.BERtarget, tx.rexdB, calc_noise_std, sim.verbose);     
end

% Transmitted power
Ptx = dBm2Watt(tx.PtxdBm);

ber.count = zeros(size(Ptx));
ber.est = zeros(size(Ptx));
ber.gauss = zeros(size(Ptx));
for k = 1:length(Ptx)
    tx.Ptx = Ptx(k);
         
    % Montecarlo simulation
    ber.count(k) = ber_apd_montecarlo(mpam, tx, fiber, apd, rx, sim);
    
    % approx ber
    [~, ber.gauss(k)] = ber_apd_doubly_stochastic(mpam, tx, fiber, apd, rx, sim);
end

if sim.verbose
    figure, hold on
    plot(tx.PtxdBm, log10(ber.count), '-o')
%     plot(tx.PtxdBm, log10(ber.est))
    plot(tx.PtxdBm, log10(ber.gauss))
    % plot(tx.PtxdBm, log10(ber.est_pdf))
    legend('Counted', 'Gaussian Approximation',...
        'Location', 'SouthWest')
    axis([tx.PtxdBm(1) tx.PtxdBm(end) -10 0])
end
    