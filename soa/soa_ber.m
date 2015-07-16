%% Calculate BER of amplified IM-DD system. 
% BER is calculated via montecarlo simulation and through saddlepoint approximation
% to calculate tail probabilities given the moment generating funcion.

function [ber, mpam] = soa_ber(mpam, tx, fiber, soa, rx, sim)

dBm2Watt = @(x) 1e-3*10.^(x/10);

% Level Spacing Optimization
if strcmp(mpam.level_spacing, 'optimized')
    % is sim.stats is set to Gaussian, then use Gaussian approximation,
    % otherwise uses accurate statistics
    if isfield(sim, 'stats') && strcmp(sim.stats, 'gaussian')
        % Optimize level spacing using Gaussian approximation
        Deltaf = rx.elefilt.noisebw(sim.fs)/2; % electric filter one-sided noise bandwidth
        Deltafopt = rx.optfilt.noisebw(sim.fs); % optical filter two-sided noise bandwidth
        % function to calculate noise std
        varTherm = rx.N0*Deltaf; % variance of thermal noise
        calc_noise_std = @(Plevel) sqrt(varTherm + 2*Plevel*soa.N0*Deltaf + 2*soa.N0^2*Deltafopt*Deltaf*(1-1/(2*Deltafopt/Deltaf)));
        % Note: Plevel corresponds to the level after SOA amplification.
        % Therefore, the soa.Gain doesn't appear in the second term because
        % it's already included in the value of Plevel.
        % Note: second term corresponds to sig-sp beat noise, and third term
        % corresponds to sp-sp beat noise with noise in one polarization.
        % Change the 2 to 4 in third term to simulate noise in two pols.
        
        % Optimize levels using Gaussian approximation
        mpam.optimize_level_spacing_gauss_approx(sim.BERtarget, tx.rexdB, calc_noise_std, sim.verbose);     
    else
        % Optimize levels using accurate noise statisitics
        [a, b] = level_spacing_optm(mpam, tx, soa, rx, sim);
        mpam.set_levels(a, b);
    end
end   
    
% Transmitted power
Ptx = dBm2Watt(tx.PtxdBm);

ber.count = zeros(size(Ptx));
ber.est = zeros(size(Ptx));
ber.gauss = zeros(size(Ptx));
for k = 1:length(Ptx)
    tx.Ptx = Ptx(k);
         
    % Montecarlo simulation
    ber.count(k) = ber_soa_montecarlo(mpam, tx, fiber, soa, rx, sim);
    
    % Estimated BER using KLSE Fourier and saddlepoint approximation of
    % tail probabilities
    [ber.est(k), ber.gauss(k)] = ber_soa_klse_fourier(mpam, tx, fiber, soa, rx, sim);
end

if sim.verbose   
    figure, hold on
    plot(tx.PtxdBm, log10(ber.count), '-o')
    plot(tx.PtxdBm, log10(ber.est))
    plot(tx.PtxdBm, log10(ber.gauss))
    % plot(tx.PtxdBm, log10(ber.est_pdf))
    legend('Counted', 'KLSE Fourier & Saddlepoint Approx', 'Gaussian Approximation',...
        'Location', 'SouthWest')
    axis([tx.PtxdBm(1) tx.PtxdBm(end) -10 0])
end
    