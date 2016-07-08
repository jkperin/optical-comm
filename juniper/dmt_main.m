%% Evaluation of DMT in pre-amplified systems
clear, close all

addpath f/
addpath ../ofdm/
addpath ../ofdm/f/
addpath ../f/
addpath ../soa/
addpath ../apd/

%% Simulation parameters
sim.Rb = 56e9;    % bit rate in bits/sec
sim.Nsymb = 2^11; % Number of symbols in montecarlo simulation
sim.Mct = 5;      % Oversampling ratio to simulate continuous time. Must be integer multiple of sim.ros.txDSP and numerator of sim.ros.rxDSP
sim.BERtarget = 1e-4; 
sim.Ndiscard = 128; % number of symbols to be discarded from the begining and end of the sequence
sim.Modulator = 'MZM'; % 'MZM' only
 
sim.RIN = true; % include RIN noise in montecarlo simulation
sim.phase_noise = true; % whether to simulate laser phase noise
sim.PMD = false; % whether to simulate PMD
sim.quantiz = ~true; % whether quantization is included
sim.stopSimWhenBERReaches0 = true; % stop simulation when counted BER reaches 0

% Control what should be plotted
sim.Plots = containers.Map();
sim.Plots('BER') = 1;
sim.Plots('Equalizer') = 0;
sim.Plots('Adaptation MSE') = 0;
sim.Plots('Constellations') = 0;
sim.Plots('Power allocation') = 0;
sim.Plots('Cyclic prefix') = 0;
sim.Plots('Estimated SNR') = 0;
sim.Plots('Channel frequency response') = 0;
sim.shouldPlot = @(x) sim.Plots.isKey(x) && sim.Plots(x);

%% OFDM 
% OFDM constructor: ofdm(Nc, Nu, CS, Rb)
% Nc : number of subcarriers
% Nu : number of used subcarriers
% CS : nominal constellation size
% Rb : bit rate (b/s)
ofdm = ofdm(512, 416, 16, 56e9, 'palloc'); 
ofdm.set_cyclic_prefix(3, 3);

%% Time and frequency
sim.N = sim.Mct*(ofdm.Nc + ofdm.Npre_os)*sim.Nsymb;           % total number of points simulated in continuous time
sim.fs = ofdm.fs*sim.Mct;
df = sim.fs/sim.N;
f = -sim.fs/2:df:sim.fs/2-df;
sim.f = f;

%% Transmitter
Tx.PtxdBm = -28:-15; % transmitter power range
% Tx.PtxdBm = -20; % transmitter power range

Tx.rexdB = -22;  % extinction ratio in dB. Defined as Pmin/Pmax
Tx.alpha = 0; % chirp parameter for laser or modulator

%% DAC
Tx.DAC.fs = ofdm.fs; % DAC sampling rate
Tx.DAC.ros = 1; % oversampling ratio of transmitter DSP
Tx.DAC.resolution = 5; % DAC effective resolution in bits
Tx.DAC.filt = design_filter('butter', 5, 30e9/(sim.fs/2)); % DAC analog frequency response
TX.DAC.rclip = 0.05;

%% Laser
% Laser constructor: laser(lambda (nm), PdBm (dBm), RIN (dB/Hz), linewidth (Hz), frequency offset (Hz))
% lambda : wavelength (nm)
% PdBm : output power (dBm)
% RIN : relative intensity noise (dB/Hz)
% linewidth : laser linewidth (Hz)
% freqOffset : frequency offset with respect to wavelength (Hz)
Tx.Laser = laser(1550e-9, 0, -150, 0.2e6, 0);

%% Modulator
Tx.Mod.type = sim.Modulator;    
Tx.Mod.BW = 25e9;
Tx.Mod.fc = Tx.Mod.BW/sqrt(sqrt(2)-1); % converts to relaxation frequency
Tx.Mod.grpdelay = 2/(2*pi*Tx.Mod.fc);  % group delay of second-order filter in seconds
Tx.Mod.Vbias = 0.5; % bias voltage normalized by Vpi
Tx.Mod.Vswing = 0.6;  % normalized voltage swing. 1 means that modulator is driven at full scale
Tx.Mod.H = @(f) exp(1j*2*pi*f*Tx.Mod.grpdelay)./(1 + 2*1j*f/Tx.Mod.fc - (f/Tx.Mod.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)

%% Fiber
% fiber(Length in m, anonymous function for attenuation versus wavelength
% (default: att(lamb) = 0 i.e., no attenuation), anonymous function for 
% dispersion versus wavelength (default: SMF28 with lamb0 = 1310nm, S0 = 0.092 s/(nm^2.km))
SMF = fiber(10e3); 
DCF = fiber(0, @(lamb) 0, @(lamb) -0.1*(lamb-1550e-9)*1e3 - 40e-6); 

Fibers = [SMF DCF];

%% Amplifier
% Class SOA characterizes amplifier in terms of gain and noise figure
% soa(GaindB: amplifier gain in dB, NFdB: noise figure in dB, lambda: operationa wavelength, maxGaindB: maximum amplifier gain)
EDFA = soa(20, 5, 1550e-9, 20); 

%% Photodiode
% pin(R: Responsivity (A/W), Id: Dark current (A), BW: bandwidth (Hz))
% PIN frequency response is modeled as a first-order system
Rx.PD = pin(1, 10e-9, Inf);

%% Receiver
% One-sided thermal noise PSD
Rx.N0 = (30e-12).^2; 

% Optical Bandpass Filter (fiber Brag gratting)
Rx.optfilt = design_filter('fbg', 4, 50e9/(sim.fs/2));

%% ADC
Rx.ADC.ros = 1;
Rx.ADC.fs = ofdm.fs;
Rx.ADC.filt = design_filter('butter', 5, 0.5*Rx.ADC.fs/(sim.fs/2)); % Antialiasing filter
Rx.ADC.ENOB = 5; % effective number of bits. Quantization is only included if sim.quantiz = true and ENOB ~= Inf
Rx.ADC.rclip = 0.05;

%% Equalizer
Rx.AdEq.mu = 1e-3;
Rx.AdEq.Ntrain = 512; % Number of symbols used in training (if Inf all symbols are used)

%% Generate summary
% generate_summary(mpam, Tx, Fibers, EDFA, Rx, sim);

%% Run simulation
ber = dmt_ber(ofdm, Tx, Fibers, EDFA, Rx, sim);