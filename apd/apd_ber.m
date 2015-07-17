%% Calculate BER of unamplified IM-DD system with APD detector 
% BER is calculated via montecarlo simulation and through saddlepoint approximation
% to calculate tail probabilities given the moment generating funcion.

function [ber, mpam] = apd_ber(mpam, tx, fiber, apd, rx, sim)

dBm2Watt = @(x) 1e-3*10.^(x/10);

% Auxiliary variables
link_gain = apd.Gain*fiber.link_attenuation(tx.lamb)*apd.R; % Overall link gain
Deltaf = rx.elefilt.noisebw(sim.fs)/2; % electric filter one-sided noise bandwidth
% function to calculate noise std
varTherm = rx.N0*Deltaf; % variance of thermal noise

if sim.RIN
    varRIN =  @(Plevel) 10^(tx.RIN/10)*Plevel.^2*Deltaf;
else
    varRIN = @(Plevel) 0;
end

% Noise std for the level Plevel
noise_std = @(Plevel) sqrt(varTherm + varRIN(Plevel) + apd.var_shot(Plevel/apd.Gain, Deltaf));

% Level Spacing Optimization
if strcmp(mpam.level_spacing, 'optimized')
    % Optimize levels using Gaussian approximation
    mpam.optimize_level_spacing_gauss_approx(sim.BERtarget, tx.rexdB, noise_std, sim.verbose);     
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
    
    % AWGN  
    mpam.adjust_levels(tx.Ptx*link_gain, tx.rexdB);

    ber.awgn(k) = mpam.ber_awgn(noise_std);
end

if sim.verbose
    figure, hold on
    plot(tx.PtxdBm, log10(ber.count), '-o')
%     plot(tx.PtxdBm, log10(ber.est))
    plot(tx.PtxdBm, log10(ber.gauss))
    % plot(tx.PtxdBm, log10(ber.est_pdf))
    plot(tx.PtxdBm, log10(ber.awgn))
    legend('Counted', 'Gaussian Approximation', 'AWGN Approximation',...
        'Location', 'SouthWest')
    axis([tx.PtxdBm(1) tx.PtxdBm(end) -10 0])
end
    