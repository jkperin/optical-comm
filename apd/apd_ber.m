function [ber, mpam, apd] = apd_ber(mpam, tx, fiber, apd, rx, sim)
%% Calculate BER of unamplified IM-DD system with APD detector 
% BER is calculated via montecarlo simulation, analytically, AWGN channel, 
% AWGN channel including noise enhancement penalty.

dBm2Watt = @(x) 1e-3*10.^(x/10);

% If equalizer is not defined assume no equalization
if ~isfield(rx, 'eq')
    rx.eq.type = 'None';
end

% Optimize APD gain
if isfield(sim, 'OptimizeGain') && sim.OptimizeGain
    apd.Gain = apd.optGain(mpam, tx, fiber, rx, sim, 'margin');
end

%% Channel response
% Hch does not include transmitter or receiver filter
if isfield(tx, 'modulator')
    Hch = tx.modulator.H(sim.f).*exp(1j*2*pi*sim.f*tx.modulator.grpdelay)...
    .*fiber.H(sim.f, tx).*apd.H(sim.f);
else
    Hch = fiber.H(sim.f, tx).*apd.H(sim.f);
end

link_gain = apd.Gain*apd.R*fiber.link_attenuation(tx.lamb); % Overall link gain

% Design equalizer
[~, eq] = equalize(rx.eq, [], Hch, mpam, rx, sim); % design equalizer
% This design assumes fixed zero-forcing equalizers

%% Noise calculations
noise_std = apd.stdNoise(eq.Hrx, eq.Hff(sim.f/mpam.Rs), rx.N0, tx.RIN, sim);

%% Level Spacing Optimization
if mpam.optimize_level_spacing  
    % Doesn't include whitening filter
    noise_std = apd.stdNoise(eq.Hrx, eq.Hff(sim.f/mpam.Rs), rx.N0, tx.RIN, sim);
    
    % Optimize levels using Gaussian approximation
    mpam = mpam.optimize_level_spacing_gauss_approx(sim.BERtarget, tx.rexdB, noise_std, sim.verbose);     
end

rx.noise_std = noise_std;

%% Calculations BER
% Transmitted power
Ptx = dBm2Watt(tx.PtxdBm);

ber.count = zeros(size(Ptx));
%ber.est = zeros(size(Ptx));
ber.gauss = zeros(size(Ptx));
ber.awgn = zeros(size(Ptx));
ber.gauss_levels = zeros(mpam.M, length(Ptx));
ber.awgn_levels = zeros(mpam.M, length(Ptx));
ber.awgn_ne = zeros(size(Ptx));
for k = 1:length(Ptx)
    tx.Ptx = Ptx(k);
            
    % Montecarlo simulation
    ber.count(k) = ber_apd_montecarlo(mpam, tx, fiber, apd, rx, sim);
    
    % Analytical BER
    [ber.gauss(k), ber.gauss_levels(:, k)] = ber_apd_gauss(mpam, tx, fiber, apd, rx, sim);
    
    % AWGN  
    mpam = mpam.adjust_levels(tx.Ptx*link_gain, tx.rexdB);

    [ber.awgn(k), ber.awgn_levels(:, k)] = mpam.ber_awgn(noise_std);
    
    % AWGN including noise enhancement penalty
    ber.awgn_ne(k) = mpam.ber_awgn(@(P) sqrt(eq.Kne)*noise_std(P));
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
    