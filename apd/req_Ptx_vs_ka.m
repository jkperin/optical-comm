%% Required Transmitted power vs impact ionization factor (ka)
function req_Ptx_vs_ka()
clc, close all

addpath ../mpam
addpath ../f
addpath f

% ka = 0.05:0.05:1;
ka = 0; % pin

% Gapd = [5:15];
Gapd = 1; % pin
M = [4 8 16];


for m = 1:length(M)
    fprintf('-- %d-PAM\n', M(m));
    for n = 1:length(Gapd)
        fprintf('---- G = %1.f\n', Gapd(n));
        [PtxdBm_req{m}(n), ber{m}(n)] = iterate(M(m), Gapd(n), false);
    end
    
%     [PtxdBm_req{m}(end+1), ber{m}(end+1), Gapd_opt{m}] = iterate(M(m), [], true);
    
%     filename = sprintf('partial_results_%dPAM.mat', M(m));
%     save(filename, 'ka', 'M', 'Gapd', 'PtxdBm_req', 'ber', 'Gapd_opt');
    1;
end

1;
% save('results_apd.mat', 'ka', 'M', 'Gapd', 'PtxdBm_req', 'ber', 'Gapd_opt');


function [PtxdBm_req, ber, Gapd_opt] = iterate(M, Gapd, OptimizeGain)
    % Simulation parameters
    sim.Nsymb = 2^10; % Number of symbols in montecarlo simulation
    sim.Mct = 9;        % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
    sim.L = 2;        % de Bruijin sub-sequence length (ISI symbol length)
    sim.BERtarget = 1e-4; 
    sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
    sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

    %
    sim.shot = true; % include shot noise. Only included in montecarlo simulation (except for APD)
    sim.RIN = true; % include RIN noise. Only included in montecarlo simulation
    sim.verbose = ~true; % show stuff

    % M-PAM
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
    switch mpam.M
        case 4
            tx.PtxdBm = -30:1:-12;
        case 8
            tx.PtxdBm = -22:1:-4;
        case 16
           tx.PtxdBm = -18:1:-2;
    end

    tx.lamb = 1310e-9; % wavelength
    tx.alpha = 2; % chirp parameter
    tx.RIN = -150;  % dB/Hz
    tx.rexdB = -10;  % extinction ratio in dB. Defined as Pmin/Pmax

    % Modulator frequency response
    % tx.modulator.fc = 2*mpam.Rs; % modulator cut off frequency
    % tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
    % tx.modulator.h = @(t) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0));
    % tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

    %% Fiber
    b2b = fiber();

    %% Receiver
    rx.N0 = (20e-12).^2; % thermal noise psd
    % Electric Lowpass Filter
    % rx.elefilt = design_filter('bessel', 5, 1.25*mpam.Rs/(sim.fs/2));
    rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);

    rx.R = 1;
    rx.Id = 10e-9;
    rx.GBw = Inf;

    %% APD 
    % (GaindB, ka, GainBW, R, Id)  
    apdG = apd(10, 1, rx.GBw, rx.R, rx.Id); 

    for k = 1:length(ka)
        apdG.ka = ka(k);

        %% Equally-spaced levels
        disp('--- Equally-spaced levels')
        mpam.level_spacing = 'equally-spaced';
        
        if OptimizeGain
            apdG.optimize_gain(mpam, tx, b2b, rx, sim);
            Gapd_opt.eq_spaced(k) = apdG.Gain;
        else
            apdG.Gain = Gapd;
            Gapd_opt = [];
        end
       
        ber_apd = apd_ber(mpam, tx, b2b, apdG, rx, sim);

        ber.eq_spaced.count(k, :) = ber_apd.count;
        ber.eq_spaced.awgn(k, :) = ber_apd.awgn;
        ber.eq_spaced.gauss(k, :) = ber_apd.gauss;  

        PtxdBm_req.eq_spaced(k) = interp1(log10(ber.eq_spaced.gauss(k, :)), tx.PtxdBm, log10(sim.BERtarget));
        
        %% Equally-spaced levels
        disp('--- Optimized Levels')
        mpam.level_spacing = 'optimized';
        
        if OptimizeGain
            apdG.optimize_gain(mpam, tx, b2b, rx, sim);
            Gapd_opt.optimized(k) = apdG.Gain;
        else
            apdG.Gain = Gapd;
            Gapd_opt = [];
        end
        
        ber_apd = apd_ber(mpam, tx, b2b, apdG, rx, sim);

        ber.optimized.count(k, :) = ber_apd.count;
        ber.optimized.awgn(k, :) = ber_apd.awgn;
        ber.optimized.gauss(k, :) = ber_apd.gauss;  

        PtxdBm_req.optimized(k) = interp1(log10(ber.optimized.gauss(k, :)), tx.PtxdBm, log10(sim.BERtarget));
    end

    figure(M), hold on, grid on
    hplot = plot(ka, PtxdBm_req.eq_spaced);
    plot(ka, PtxdBm_req.optimized, '--', 'Color', get(hplot, 'Color'))
    xlabel('k_a')
    ylabel('Required Transmitted Optical Power (dBm)')
    legend('Equally-spaced', 'Optimized')
    % axis([0 1 -21 -19]);

    if OptimizeGain
        figure(M+1), hold on, grid on
        hplot = plot(ka, Gapd_opt.eq_spaced);
        plot(ka, Gapd_opt.optimized, '--', 'Color', get(hplot, 'Color'))
        xlabel('k_a')
        ylabel('Optimal APD Gain')
        legend('Equally-spaced', 'Optimized')
    end
end
end