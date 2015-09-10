function [ber, mpam] = soa_ber(mpam, tx, fiber, soa, rx, sim)
%% Calculate BER of amplified IM-DD system with SOA. 
% BER is calculated via montecarlo simulation, analytically using saddlepoint
% approximation, AWGN channel, AWGN channel including noise enhancement penalty.

dBm2Watt = @(x) 1e-3*10.^(x/10);

% Channel response
Ptx = design_filter('matched', mpam.pshape, 1/sim.Mct); % transmitted pulse shape

% Hch does not include receiver filter
if isfield(tx, 'modulator')
    Hch = Ptx.H(sim.f/sim.fs).*tx.modulator.H(sim.f).*exp(1j*2*pi*sim.f*tx.modulator.grpdelay)...
    .*fiber.H(sim.f, tx);
else
    Hch = Ptx.H(sim.f/sim.fs).*fiber.H(sim.f, tx);
end

link_gain = soa.Gain*rx.R*fiber.link_attenuation(tx.lamb); % Overall link gain

%% Polarizer
if isfield(sim, 'polarizer') && ~sim.polarizer
    Npol = 2;     % number of noise polarizations
else % by default assumes that polarizer is being use so Npol = 1.
    Npol = 1;
end

%% Noise calculations
% Thermal noise
Deltaf = rx.elefilt.noisebw(sim.fs)/2; % electric filter one-sided noise bandwidth
Deltafopt = rx.optfilt.noisebw(sim.fs); % optical filter noise bandwidth
% Thermal noise
varTherm = rx.N0*Deltaf; % variance of thermal noise

% Shot noise
if isfield(sim, 'shot') && sim.shot
    varShot = @(Plevel) 2*1.60217657e-19*(rx.R*Plevel + rx.Id)*Deltaf;
else
    varShot = @(Plevel) 0;
end

% RIN
if isfield(sim, 'RIN') && sim.RIN
    varRIN =  @(Plevel) 10^(tx.RIN/10)*Plevel.^2*Deltaf;
else
    varRIN = @(Plevel) 0;
end

% Noise std for the level Plevel
noise_std = @(Plevel) sqrt(varTherm + varShot(Plevel) + rx.R^2*varRIN(Plevel)...
    + rx.R^2*soa.var_awgn(Plevel/soa.Gain, Deltaf, Deltafopt, Npol));
% Note: Plevel is divided by SOA gain to obtain power at the amplifier input

% Noise enhancement penalty
if isfield(rx, 'eq')
    [~, rx.eq] = equalize(rx.eq, [], Hch, mpam, rx, sim); % design equalizer
    % This design assumes fixed zero-forcing equalizers
    Kne = rx.eq.Kne; % noise enhancement penalty
    % Kne = noise variance after equalizer/noise variance before equalizer
else 
    rx.eq.type = 'None';
    Kne = 1;
end

%% Level Spacing Optimization
if strcmp(mpam.level_spacing, 'optimized')
    % is sim.stats is set to Gaussian, then use Gaussian approximation,
    % otherwise uses accurate statistics
    if isfield(sim, 'stats') && strcmp(sim.stats, 'gaussian')        
        % Optimize levels using Gaussian approximation
        mpam = mpam.optimize_level_spacing_gauss_approx(sim.BERtarget, tx.rexdB, noise_std, sim.verbose);     
    else
        % Optimize levels using accurate noise statisitics
        [a, b] = level_spacing_optm(mpam, tx, soa, rx, sim);
        mpam = mpam.set_levels(a, b);
    end
end   

%% Calculations BER
% Transmitted power
Ptx = dBm2Watt(tx.PtxdBm);

ber.count = zeros(size(Ptx));
ber.est = zeros(size(Ptx));
ber.gauss = zeros(size(Ptx));
ber.awgn = zeros(size(Ptx));
ber.awgn_ne = zeros(size(Ptx));
for k = 1:length(Ptx)
    tx.Ptx = Ptx(k);
         
    % Montecarlo simulation
    ber.count(k) = ber_soa_montecarlo(mpam, tx, fiber, soa, rx, sim);
    
    % Estimated BER using KLSE Fourier and saddlepoint approximation of
    % tail probabilities
    [ber.est(k), ber.gauss(k)] = ber_soa_klse_fourier(mpam, tx, fiber, soa, rx, sim);
    
    % AWGN  
    mpam = mpam.adjust_levels(tx.Ptx*link_gain, tx.rexdB);

    ber.awgn(k) = mpam.ber_awgn(noise_std);
    
    % AWGN including noise enhancement penalty
    ber.awgn_ne(k) = mpam.ber_awgn(@(P) sqrt(Kne)*noise_std(P));
end

if sim.verbose   
    figure, hold on
    plot(tx.PtxdBm, log10(ber.count), '-o')
    plot(tx.PtxdBm, log10(ber.est))
    plot(tx.PtxdBm, log10(ber.gauss))
    plot(tx.PtxdBm, log10(ber.awgn))
    plot(tx.PtxdBm, log10(ber.awgn_ne))
    % plot(tx.PtxdBm, log10(ber.est_pdf))
    legend('Counted', 'KLSE Fourier & Saddlepoint Approx', 'Gaussian approximation', 'AWGN',...
        'AWGN + noise enhancement', 'Location', 'Location', 'SouthWest')
    axis([tx.PtxdBm(1) tx.PtxdBm(end) -10 0])
end
    