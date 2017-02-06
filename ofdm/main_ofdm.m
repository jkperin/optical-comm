%% Evaluation of OFDM in IM-DD system, which may be amplified or not
clear, close all

addpath f/
addpath ../f/
addpath ../apd/

%% Transmit power swipe
Tx.PtxdBm = -30:-22; % transmitter power range
Tx.PtxdBm = -10:-5; % transmitter power range

%% Simulation parameters
sim.Rb = 56e9;    % bit rate in bits/sec
sim.Nsymb = 2^10; % Number of symbols in montecarlo simulation
sim.Mct = 5;      % Oversampling ratio to simulate continuous time. Must be integer multiple of sim.ros.txDSP and numerator of sim.ros.rxDSP
sim.BERtarget = 1e-4; 
sim.Ndiscard = 32; % number of symbols to be discarded from the begining and end of the sequence
sim.Modulator = 'DML'; % 'MZM' or 'DML'
 
%% Simulation control
sim.preAmp = false;
sim.preemphasis = false; % preemphasis to compensate for transmitter bandwidth limitation
sim.preemphRange = 25e9; % preemphasis range
sim.RIN = true; % include RIN noise in montecarlo simulation
sim.phase_noise = true; % whether to simulate laser phase noise
sim.PMD = false; % whether to simulate PMD
sim.quantiz = true; % whether quantization is included
sim.stopSimWhenBERReaches0 = true; % stop simulation when counted BER reaches 0

% Control what should be plotted
sim.Plots = containers.Map();
sim.Plots('BER') = 1;
sim.Plots('Equalizer') = 1;
sim.Plots('Adaptation MSE') = 0;
sim.Plots('Constellations') = 1;
sim.Plots('Power allocation') = 1;
sim.Plots('Cyclic prefix') = 1;
sim.Plots('Estimated SNR') = 1;
sim.Plots('Channel frequency response') = 0;
sim.Plots('OSNR') = 1;
sim.shouldPlot = @(x) sim.Plots.isKey(x) && sim.Plots(x);

%% OFDM 
% OFDM constructor: ofdm(Nc, Nu, CS, Rb, power_allocation_type (optional, default = 'palloc')
% Nc : number of subcarriers
% Nu : number of used subcarriers
% CS : nominal constellation size
% Rb : bit rate (b/s)
% power_allocation_type : {'palloc', 'preemphasis'}
ofdm = ofdm(512, 416, 16, 56e9, 'preemphasis'); 
ofdm.set_cyclic_prefix(3, 3); % set cyclic prefix length. Should be consistent with channel memory length

%% Time and frequency
sim.N = sim.Mct*(ofdm.Nc + ofdm.Npre_os)*sim.Nsymb;           % total number of points simulated in continuous time
sim.fs = ofdm.fs*sim.Mct;
[sim.f, sim.t] = freq_time(sim.N, sim.fs);

%% ==========================  Transmitter ================================
%% DAC
Tx.DAC.fs = ofdm.fs; % DAC sampling rate
Tx.DAC.ros = 1; % oversampling ratio of transmitter DSP
Tx.DAC.resolution = 5; % DAC effective resolution in bits
Tx.DAC.filt = design_filter('butter', 5, 25e9/(sim.fs/2)); % DAC analog frequency response
% TX.DAC.rclip = 0.05;
% juniperDACfit = [-0.0013 0.5846 0];
% Tx.DAC.filt.H = @(f) 1./(10.^(polyval(juniperDACfit, abs(f*sim.fs/1e9))/20)); % DAC analog frequency response

%% Modulator
Tx.rexdB = -20;  % extinction ratio in dB. Defined as Pmin/Pmax
Tx.rclip = 4; % clipping ratio
Tx.alpha = 0; % chirp parameter
Tx.Mod.type = sim.Modulator;    
Tx.Mod.BW = 25e9;
Tx.Mod.filt = design_filter('bessel', 5, Tx.Mod.BW/(sim.fs/2));
Tx.Mod.H = Tx.Mod.filt.H(sim.f/sim.fs);

%% Laser
% Laser constructor: laser(lambda (nm), PdBm (dBm), RIN (dB/Hz), linewidth (Hz), frequency offset (Hz))
% lambda : wavelength (nm)
% PdBm : output power (dBm)
% RIN : relative intensity noise (dB/Hz)
% linewidth : laser linewidth (Hz)
% freqOffset : frequency offset with respect to wavelength (Hz)
Tx.Laser = laser(1550e-9, 0, -150, 0.2e6, 0);
Tx.Laser.H = @(f) Tx.Mod.filt.H(f/sim.fs); % only used if Tx.Mod.type = 'DML'

%% Fiber
% fiber(Length in m, anonymous function for attenuation versus wavelength
% (default: att(lamb) = 0 i.e., no attenuation), anonymous function for 
% dispersion versus wavelength (default: SMF28 with lamb0 = 1310nm, S0 = 0.092 s/(nm^2.km))
SMF = fiber(0e3); 
DCF = fiber(0, @(lamb) 0, @(lamb) -0.1*(lamb-1550e-9)*1e3 - 40e-6); 

Fibers = [SMF DCF];

linkAttdB = SMF.att(Tx.Laser.wavelength)*SMF.L/1e3...
    + DCF.att(Tx.Laser.wavelength)*DCF.L/1e3;

%% ========================== Amplifier ===================================
% Constructor: OpticalAmplifier(GaindB, NF, lamb, maxGain (optional))
% GaindB : Gain in dB
% NF : Noise figure in dB
% lamb : wavelength in m
% maxGain = maximum amplifier gain in dB, default is Inf
Rx.OptAmp = OpticalAmplifier(30, 5, Tx.Laser.lambda);
Rx.OptAmpOutPowerdBm = 0; % output power after amplifier
% Note: the amplifier here operates in the constant output power mode,
% where the output power after amplification is set to Rx.AmpOutPowerdBm

%% ============================ Receiver ==================================
%% Photodiode
% pin(R: Responsivity (A/W), Id: Dark current (A), BW: bandwidth (Hz))
% PIN frequency response is modeled as a first-order system
Rx.PD = pin(1, 10e-9, Inf);

%% TIA-AGC
% One-sided thermal noise PSD
Rx.N0 = (30e-12).^2; 

%% Receiver DSP
Rx.filtering = 'antialiasing'; % {'antialiasing' or 'matched'}

%% ADC
Rx.ADC.ros = 1;
Rx.ADC.fs = ofdm.fs;
Rx.ADC.filt = design_filter('butter', 5, 0.5*Rx.ADC.fs/(sim.fs/2)); % Antialiasing filter
Rx.ADC.ENOB = 5; % effective number of bits. Quantization is only included if sim.quantiz = true and ENOB ~= Inf
% Rx.ADC.rclip = 0.05;

%% Equalizer
Rx.AdEq.mu = 1e-2;
Rx.AdEq.Ntrain = 512; % Number of frames used in training (if Inf all symbols are used)

%% Generate summary
% generate_summary(mpam, Tx, Fibers, EDFA, Rx, sim);

%% Run simulation
[berOFDM, ofdm, SNRdB, OSNRdB] = ber_ofdm(ofdm, Tx, Fibers, Rx, sim)