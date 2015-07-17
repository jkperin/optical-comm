%% Required Transmitted power vs impact ionization factor (ka)
function req_Ptx_vs_ka()
clc, close all

addpath ../mpam
addpath ../f
addpath f

ka = 0.01:0.05:1;

Gapd = [1:15];
M = [4 8 16];


for m = 1:length(M)
    for n = 1:length(Gapd)
        [PtxdBm_req{m}(n, :), ber{m}(n)] = iterate(M(m), Gapd(n));
    end
    
     [PtxdBm_req{m}(end+1, :), ber{m}(end+1), Gapd_opt{m}] = iterate(M(m), []);
end
1;

function [PtxdBm_req, ber, Gapd_opt] = iterate(M, Gapd)
    % Simulation parameters
    sim.Nsymb = 2^15; % Number of symbols in montecarlo simulation
    sim.Mct = 9;        % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
    sim.L = 2;        % de Bruijin sub-sequence length (ISI symbol length)
    sim.BERtarget = 1e-4; 
    sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
    sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

    %
    sim.shot = true; % include shot noise. Only included in montecarlo simulation (except for APD)
    sim.RIN = ~true; % include RIN noise. Only included in montecarlo simulation
    sim.verbose = ~true; % show stuff

    % M-PAM
    mpam = PAM(M, 100e9, 'optimized', @(n) double(n >= 0 & n < sim.Mct));

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
            tx.PtxdBm = -22:2:-4;
        case 16
           tx.PtxdBm = -18:2:-2;
    end

    tx.lamb = 1310e-9; % wavelength
    tx.alpha = 2; % chirp parameter
    tx.RIN = -150;  % dB/Hz
    tx.rexdB = -Inf;  % extinction ratio in dB. Defined as Pmin/Pmax

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

    ber.count = zeros(length(ka), length(tx.PtxdBm));
    ber.awgn = zeros(length(ka), length(tx.PtxdBm));
    ber.gauss = zeros(length(ka), length(tx.PtxdBm));
    Gapd_opt = zeros(size(ka));
    PtxdBm_req = zeros(size(ka));
    for k = 1:length(ka)
        apdG.ka = ka(k);

        if isempty(Gapd)
            apdG.optimize_gain(mpam, tx, b2b, rx, sim);
            Gapd_opt(k) = apdG.Gain;
        else
            apdG.Gain = Gapd;
            Gapd_opt = [];
        end

        ber_apd = apd_ber(mpam, tx, b2b, apdG, rx, sim);

        ber.count(k, :) = ber_apd.count; 
        ber.awgn(k, :) = ber_apd.awgn;
        ber.gauss(k, :) = ber_apd.gauss;  

        PtxdBm_req(k) = interp1(log10(ber_apd.gauss), tx.PtxdBm, log10(sim.BERtarget));

        % 
    end

    figure(1), hold on, grid on
    plot(ka, PtxdBm_req)
    xlabel('k_a')
    ylabel('Required Transmitted Optical Power (dBm)')
    % axis([0 1 -21 -19]);

    if ~isempty(Gapd_opt)
        figure(2), hold on, grid on
        plot(ka, Gapd_opt)
        xlabel('k_a')
        ylabel('Optimal APD Gain')
    end
end
end