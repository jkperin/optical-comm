%% Required Transmitted power vs amplifier noise figure
function req_Ptx_vs_Fn()
clc, close all

addpath ../mpam/
addpath ../f % general functions
addpath f

%% Parameters to swipe
M = [4 8 16]; % PAM order
Fn = 3:10; % Noise figure
GainsdB = 20; % Gain in dB

for m = 1:length(M)
    fprintf('------- %d-PAM\n (with polarizer)', M(m))
    
    [PtxdBm_req_polarizer{m}, ber_polarizer{m}] = iterate(M(m), Fn, GainsdB, true);
    
    fprintf('------- %d-PAM\n (no polarizer)', M(m))
    [PtxdBm_req{m}, ber{m}] = iterate(M(m), Fn, GainsdB, false);
   
end

save soa_awgn_results_req_Ptx_vs_Fn M Fn GainsdB PtxdBm_req PtxdBm_req_polarizer ber ber_polarizer

function [PtxdBm_req, ber] = iterate(M, Fn, GaindB, polarizer)

    % Simulation parameters
    sim.Nsymb = 2^10; % Number of symbols in montecarlo simulation
    sim.Mct = 15;     % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
    sim.L = 2;        % de Bruijin sub-sequence length (ISI symbol length)
    sim.Me = 16; % Number of used eigenvalues
    sim.BERtarget = 1e-4; 
    sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
    sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

    sim.polarizer = polarizer; % if true uses 1 noise polarization
    sim.shot = true; % include shot noise in montecarlo simulation 
    sim.RIN = true; % include RIN noise in montecarlo simulation
    sim.verbose = false; % show stuff

    % M-PAM
    % M, Rb, leve_spacing, pshape
    mpam = PAM(M, 100e9, 'equally-spaced', @(n) double(n >= 0 & n < sim.Mct));

    %% Time and frequency
    sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'

    dt = 1/sim.fs;
    t = (0:dt:(sim.N-1)*dt).';
    df = 1/(dt*sim.N);
    f = (-sim.fs/2:df:sim.fs/2-df).';

    sim.t = t;
    sim.f = f;

    %% Transmitter
    if mpam.M == 8
        tx.PtxdBm = -22:2:-8;
    elseif mpam.M == 4
        tx.PtxdBm = -26:1:-10;
    elseif mpam.M == 16
        tx.PtxdBm = -18:2:-2;
    end

    tx.lamb = 1310e-9; % wavelength
    tx.alpha = 0; % chirp parameter
    tx.RIN = -150;  % dB/Hz
    tx.rexdB = -10;  % extinction ratio in dB. Defined as Pmin/Pmax

    % Modulator frequency response
    % tx.modulator.BW = 2*mpam.Rs; % modulator cut off frequency
    % tx.modulator.fc = tx.modulator.BW/sqrt(sqrt(2)-1); % converts to relaxation frequency
    % tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
    % tx.modulator.h = @(t) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0));
    % tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

    %% Receiver
    rx.N0 = (20e-12).^2; % thermal noise psd
    rx.Id = 10e-9; % dark current
    rx.R = 1; % responsivity
    % Electric Lowpass Filter
    % rx.elefilt = design_filter('bessel', 5, 0.5*mpam.Rs/(sim.fs/2));
    rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);
    % Optical Bandpass Filter
    rx.optfilt = design_filter('fbg', 0, 200e9/(sim.fs/2));

    % KLSE Fourier Series Expansion (done here because depends only on filters
    % frequency response)
    [rx.U_fourier, rx.D_fourier, rx.Fmax_fourier] = klse_fourier(rx, sim, sim.Mct*(mpam.M^sim.L + 2*sim.L)); 

    %% SOA
    % soa(GaindB, NF, lambda, maxGaindB)
    soaG = soa(GaindB, 9, tx.lamb, Inf);

    for k= 1:length(Fn)
        fprintf('---- Fn = %.1f dB\n', Fn(k));
        soaG.Fn = Fn(k); % set Noise Figure

        %% Equally-spaced levels
        disp('-- Equally-spaced levels')
        mpam.level_spacing = 'equally-spaced';

        ber_soa = soa_ber(mpam, tx, fiber, soaG, rx, sim);
        
        ber.eq_spaced.est(k, :) = ber_soa.est;  
        ber.eq_spaced.count(k, :) = ber_soa.count;
        ber.eq_spaced.awgn(k, :) = ber_soa.awgn;

        PtxdBm_req.eq_spaced(k) = interp1(log10(ber.eq_spaced.est(k, :)), tx.PtxdBm, log10(sim.BERtarget), 'spline');

        %% Optimized levels
        disp('-- Optimized levels')
        mpam.level_spacing = 'optimized';

        ber_soa = soa_ber(mpam, tx, fiber, soaG, rx, sim);
        
        ber.optimized.est(k, :) = ber_soa.est; 
        ber.optimized.count(k, :) = ber_soa.count;
        ber.optimized.awgn(k, :) = ber_soa.awgn;

        PtxdBm_req.optimized(k) = interp1(log10(ber.optimized.est(k, :)), tx.PtxdBm, log10(sim.BERtarget), 'spline');
    end

    %% Plot
    legends = {};
    for k = 1:length(Fn)
        figure(M), hold on, grid on, box on
        hplot(k) = plot(tx.PtxdBm, log10(ber.eq_spaced.count(k, :)), '-o');
        legends = [legends, sprintf('Fn = %.1f dB', Fn(k))];

        figure(M+1), hold on, grid on, box on
        plot(tx.PtxdBm, log10(ber.optimized.count(k, :)), '-o', 'Color', get(hplot(k), 'Color'));
    end

    for k = 1:length(Fn)
        figure(M)
        plot(tx.PtxdBm, log10(ber.eq_spaced.est(k, :)), '-', 'Color', get(hplot(k), 'Color'));
        plot(tx.PtxdBm, log10(ber.eq_spaced.awgn(k, :)), '--', 'Color', get(hplot(k), 'Color'));

        figure(M+1)
        plot(tx.PtxdBm, log10(ber.optimized.est(k, :)), '-', 'Color', get(hplot(k), 'Color'));
        plot(tx.PtxdBm, log10(ber.optimized.awgn(k, :)), '--', 'Color', get(hplot(k), 'Color'));
    end

    figure(M+1)
    xlabel('Received Power (dBm)')
    ylabel('log(BER)')
    legend(legends{:})
    axis([tx.PtxdBm([1 end]) -8 0])

    figure(M+1)
    xlabel('Received Power (dBm)')
    ylabel('log(BER)')
    legend(legends{:})
    axis([tx.PtxdBm([1 end]) -8 0])

    figure(1), box on, grid on, hold on
    hplot = plot(Fn, PtxdBm_req.eq_spaced, '-o');
    plot(Fn, PtxdBm_req.optimized, '-s', 'Color', get(hplot, 'Color'));
    xlabel('Noise Figure (dB)')
    ylabel('Required Transmitted Power (dBm)')
end
end