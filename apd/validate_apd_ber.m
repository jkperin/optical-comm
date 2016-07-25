%% Validate APD BER calculation using Montecarlo simulation, enummeration, and AWGN approximation including noise enhancement
clear, clc

addpath ../mpam
addpath ../f
addpath f

%% Simulation parameters
sim.Nsymb = 2^15; % Number of symbols in montecarlo simulation
sim.ros.txDSP = 1; % oversampling ratio of transmitter DSP
sim.ros.rxDSP = 1; % oversampling ratio of receivre DSP
sim.Mct = 12*sim.ros.rxDSP;    % Oversampling ratio to simulate continuous time (must be odd so that sampling is done right, and FIR filters have interger grpdelay)  
sim.L = 4;        % de Bruijin sub-sequence length (ISI symbol length)
sim.BERtarget = 1.8e-4; 
sim.Ndiscard = 256;  % number of 0 symbols to be inserted at the begining and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation
sim.WhiteningFilter = true;

%
sim.RIN = true; % include RIN noise. Only included in montecarlo simulation
sim.quantiz = false; % whether quantization is included
sim.terminateWhenBERReaches0 = true; % whether simulation is terminated when counted BER reaches 0

sim.Plots = containers.Map();
sim.Plots('BER') = 1;
% sim.Plots('Adaptation MSE') = 1;
% sim.Plots('Frequency Response') = 1;
% sim.Plots('Equalizer') = 1;
% sim.Plots('DAC output') = 1;
% sim.Plots('Optical eye diagram') = 1;
% sim.Plots('Received signal eye diagram') = 1;
% sim.Plots('Signal after equalization') = 1;
sim.shouldPlot = @(x) sim.Plots.isKey(x) && sim.Plots(x);

%% M-PAM
pulse_shape = select_pulse_shape('rect', sim.ros.txDSP);
mpam = PAM(4, 107e9, 'equally-spaced', pulse_shape);

%% Time and frequency
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'
[sim.f, sim.t] = freq_time(sim.N, sim.fs);

%% Laser
% Laser constructor: laser(lambda (nm), PdBm (dBm), RIN (dB/Hz), linewidth (Hz), frequency offset (Hz))
% lambda : wavelength (nm)
% PdBm : output power (dBm)
% RIN : relative intensity noise (dB/Hz)
% linewidth : laser linewidth (Hz)
% freqOffset : frequency offset with respect to wavelength (Hz)
Tx.Laser = laser(1250e-9, 0, -150, 0.2e6, 0);

%% Transmitter
Tx.PtxdBm = -22:-10; % transmitted power swipe
% Tx.PtxdBm = -15; % transmitted power swipe

%% DAC
Tx.DAC.fs = sim.ros.txDSP*mpam.Rs; % DAC sampling rate
Tx.DAC.ros = sim.ros.txDSP; % oversampling ratio of transmitter DSP
Tx.DAC.resolution = Inf; % DAC effective resolution in bits
Tx.DAC.filt = design_filter('bessel', 5, 200e9/(sim.fs/2)); % DAC analog frequency response

%% Modulator
Tx.Mod.alpha = 2; % chirp parameter
Tx.Mod.rexdB = -10;  % extinction ratio in dB. Defined as Pmin/Pmax
% Modulator frequency response
Tx.Mod.BW = 30e9;
Tx.Mod.fc = Tx.Mod.BW/sqrt(sqrt(2)-1); % converts to relaxation frequency
Tx.Mod.h = @(t) (2*pi*Tx.Mod.fc)^2*t(t >= 0).*exp(-2*pi*Tx.Mod.fc*t(t >= 0));
Tx.Mod.grpdelay = 2/(2*pi*Tx.Mod.fc);  % group delay of second-order filter in seconds
Tx.Mod.H = @(f) exp(1j*2*pi*f.*Tx.Mod.grpdelay)./(1 + 2*1j*f/Tx.Mod.fc - (f/Tx.Mod.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)

%% Fiber
Fiber = fiber(0e3);

%% Receiver
Rx.N0 = (30e-12).^2; % thermal noise psd
Rx.filtering = 'matched'; % Electric Lowpass Filter prior to sampling and equalization: either 'antialisaing' or 'matched'

%% ADC
Rx.ADC.ros = sim.ros.rxDSP;
Rx.ADC.fs = Rx.ADC.ros*mpam.Rs*(mpam.pulse_shape.rolloff + 1)/2;
Rx.ADC.filt = design_filter('butter', 5, 0.5*Rx.ADC.fs/(sim.fs/2)); % Antialiasing filter
Rx.ADC.ENOB = Inf; % effective number of bits. Quantization is only included if sim.quantiz = true and ENOB ~= Inf
Rx.ADC.rclip = 0;

%% Equalization
% Rx.eq.type = 'None';
% Rx.eq.type = 'Adaptive TD-LE';
Rx.eq.type = 'Fixed TD-SR-LE';
Rx.eq.ros = sim.ros.rxDSP;
Rx.eq.Ntaps = 31;
% Rx.eq.mu = 1e-3;
% Rx.eq.Ntrain = Inf; % Number of symbols used in training (if Inf all symbols are used)
Rx.eq.Ndiscard = [1024 1024]; % symbols to be discard from begining and end of sequence due to adaptation, filter length, etc

%% APD 
% apd(GaindB, ka, [BW GBP=Inf], R, Id) 
Apd = apd(10, 0.1, [20e9 270e9], 0.74, 10e-9);
% Apd = apd(15, 0.5, Inf, 1, 10e-9);

% Apd = pin(1, 1e-9);

% BER
sim.OptimizeGain = true;
ber_apd = apd_ber(mpam, Tx, Fiber, Apd, Rx, sim);
% 
mpam.level_spacing = 'optimized';
ber_apd_eq = apd_ber(mpam, Tx, Fiber, Apd, Rx, sim);
        