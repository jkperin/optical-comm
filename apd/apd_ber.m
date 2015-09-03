function [ber, mpam] = apd_ber(mpam, tx, fiber, apd, rx, sim)
%% Calculate BER of unamplified IM-DD system with APD detector 
% BER is calculated via montecarlo simulation, analytically, AWGN channel, 
% AWGN channel including noise enhancement penalty.

dBm2Watt = @(x) 1e-3*10.^(x/10);

% Auxiliary variables
link_gain = apd.Gain*fiber.link_attenuation(tx.lamb)*apd.R; % Overall link gain

%% Noise calculations
% Thermal noise
Deltaf = rx.elefilt.noisebw(sim.fs)/2; % electric filter one-sided noise bandwidth
varTherm = rx.N0*Deltaf; % variance of thermal noise

% RIN
if isfield(sim, 'RIN') && sim.RIN
    varRIN =  @(Plevel) 10^(tx.RIN/10)*Plevel.^2*Deltaf;
else
    varRIN = @(Plevel) 0;
end

% Noise std for the level Plevel
noise_std = @(Plevel) sqrt(varTherm + varRIN(Plevel) + apd.varShot(Plevel/apd.Gain, Deltaf));
% Note: Plevel is divided by APD gain to obtain power at the apd input

% Noise enhancement penalty
if isfield(rx, 'eq')
    [~, rx.eq] = equalize(rx.eq.type, [], mpam, tx, fiber, rx, sim); % design equalizer
    % This design assumes fixed zero-forcing equalizers
    Kne = rx.eq.Kne; % noise enhancement penalty
    % Kne = noise variance after equalizer/noise variance before equalizer
else 
    rx.eq.type = 'None';
    Kne = 1;
end

%% Level Spacing Optimization
if mpam.optimize_level_spacing
    % Optimize levels using Gaussian approximation
    mpam.optimize_level_spacing_gauss_approx(sim.BERtarget, tx.rexdB, noise_std, sim.verbose);     
end

%% Calculations BER
% Transmitted power
Ptx = dBm2Watt(tx.PtxdBm);

ber.count = zeros(size(Ptx));
ber.est = zeros(size(Ptx));
ber.gauss = zeros(size(Ptx));
ber.awgn = zeros(size(Ptx));
ber.awgn_levels = zeros(mpam.M, length(Ptx));
ber.awgn_ne = zeros(size(Ptx));
for k = 1:length(Ptx)
    tx.Ptx = Ptx(k);
            
    % Montecarlo simulation
    ber.count(k) = ber_apd_montecarlo(mpam, tx, fiber, apd, rx, sim);
    
    % Analytical BER
    [~, ber.gauss(k)] = ber_apd_doubly_stochastic(mpam, tx, fiber, apd, rx, sim);
    
    % AWGN  
    mpam.adjust_levels(tx.Ptx*link_gain, tx.rexdB);

    [ber.awgn(k), ber.awgn_levels(:, k)] = mpam.ber_awgn(noise_std);
    
    % AWGN including noise enhancement penalty
    ber.awgn_ne(k) = mpam.ber_awgn(@(P) sqrt(Kne)*noise_std(P));
end

if sim.verbose
    figure, hold on
    plot(tx.PtxdBm, log10(ber.count), '-o')
%     plot(tx.PtxdBm, log10(ber.est))
    plot(tx.PtxdBm, log10(ber.gauss))
    % plot(tx.PtxdBm, log10(ber.est_pdf))
    plot(tx.PtxdBm, log10(ber.awgn))
    plot(tx.PtxdBm, log10(ber.awgn_ne))
    legend('Counted', 'Gaussian approximation', 'AWGN approximation',...
        'AWGN + noise enhancement', 'Location', 'SouthWest')
    axis([tx.PtxdBm(1) tx.PtxdBm(end) -10 0])
end
    