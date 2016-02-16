function [ber, mpam, apd] = apd_ber(mpam, tx, fiber, apd, rx, sim)
%% Calculate BER of unamplified IM-DD system with APD detector 
% BER is calculated via montecarlo simulation, analytically, AWGN channel, 
% AWGN channel including noise enhancement penalty.
verbose = sim.verbose;
sim.verbose = max(sim.verbose-1, 0);

dBm2Watt = @(x) 1e-3*10.^(x/10);

% If equalizer is not defined assume no equalization
if ~isfield(rx, 'eq')
    rx.eq.type = 'None';
end

% Optimize APD gain
if isfield(sim, 'OptimizeGain') && sim.OptimizeGain
    [apd.Gain, mpam] = apd.optGain(mpam, tx, fiber, rx, sim);
    fprintf('Optimal APD Gain = %.2f (%2.f dB)\n', apd.Gain, apd.GaindB);
    % if mpam.level_spacing = 'optimized', then apd.optGain returns mpam 
    % with optimal level spacing
elseif mpam.optimize_level_spacing  %% Level Spacing Optimization
    % Optimize levels using Gaussian approximation
    [~, mpam] = apd.optimize_PAM_levels(apd.Gain, mpam, tx, fiber, rx, sim);
    mpam = mpam.norm_levels();
end

%% BER
% Transmitted power
Ptx = dBm2Watt(tx.PtxdBm);

ber.count = zeros(size(Ptx)); % counted BER
ber.gauss = zeros(size(Ptx)); % analysis assuming Gaussian stats
ber.awgn = zeros(size(Ptx)); % AWGN approximation (includes noise enhancement penalty)
ber.gauss_levels = zeros(mpam.M, length(Ptx)); % symbol error prob for each level
for k = 1:length(Ptx)
    tx.Ptx = Ptx(k);
            
    % Montecarlo simulation
    ber.count(k) = ber_apd_montecarlo(mpam, tx, fiber, apd, rx, sim);
    
    % BER using Gaussian stats approximation for shot noise (enumeration)
    % AWGN approximation is also included by this function
    [ber.gauss(k), ber.gauss_levels(:, k), ber.awgn(k)] = ...
        ber_apd_gauss(mpam, tx, fiber, apd, rx, sim);
end

if verbose
    PrxdBm = tx.PtxdBm - 10*log10(fiber.link_attenuation(tx.lamb));
    figure(1), hold on, box on
    hline = plot(PrxdBm, log10(ber.count), 'o');
    plot(PrxdBm, log10(ber.gauss), '-', 'Color', get(hline, 'Color'))
    plot(PrxdBm, log10(ber.awgn), '--', 'Color', get(hline, 'Color'))
    legend('Counted', 'Gaussian stats approximation', 'AWGN approximation',...
        'Location', 'SouthWest')
    xlabel('Received Power (dBm)')
    ylabel('log_{10}(BER)')
    grid on
    axis([tx.PtxdBm(1) tx.PtxdBm(end) -8 0])
end