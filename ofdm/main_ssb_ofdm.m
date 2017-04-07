%% Evaluation of OFDM in IM-DD system, which may be amplified or not
clear, clc, close all

addpath f/
addpath ../f/
addpath ../apd/

%% Transmit power swipe
% Tx.PtxdBm = -10:-5; % transmitter power range
% Tx.PtxdBm = -12:-7; % transmitter power range
Tx.PtxdBm = -10; % transmitter power range

%% Simulation parameters
sim.Rb = 112e9;    % bit rate in bits/sec
sim.Nsymb = 2^11; % Number of symbols in montecarlo simulation
sim.Mct = 3;      % Oversampling ratio to simulate continuous time. Must be integer multiple of sim.ros.txDSP and numerator of sim.ros.rxDSP
sim.BERtarget = 1.8e-4; 
sim.Ndiscard = 128; % number of symbols to be discarded from the begining and end of the sequence
sim.OFDM = 'SSB-OFDM'; 
 
%% Simulation control
sim.RIN = true; % include RIN noise in montecarlo simulation
sim.phase_noise = false; % whether to simulate laser phase noise
sim.PMD = false; % whether to simulate PMD
sim.quantiz = true; % whether quantization is included
sim.SSBIcancellation = true;
sim.stopSimWhenBERReaches0 = true; % stop simulation when counted BER reaches 0

% Control what should be plotted
sim.Plots = containers.Map();
sim.Plots('BER') = 1;
sim.Plots('Equalizer') = 1;
sim.Plots('Adaptation MSE') = 0;
sim.Plots('Constellations') = 1;
sim.Plots('Power allocation') = 1;
sim.Plots('Cyclic prefix') = 0;
sim.Plots('Estimated SNR') = 0;
sim.Plots('Decision errors') = 1;
sim.Plots('Channel frequency response') = 0;
sim.Plots('OSNR') = 1;
sim.shouldPlot = @(x) sim.Plots.isKey(x) && sim.Plots(x);

%% OFDM 
% OFDM constructor: ofdm(Nc, Nu, CS, Rb, prefix (optional, default = 'DSB'), power_allocation_type (optional, default = 'Levin-Campello-MA')
% Nc : number of subcarriers
% Nu : number of used subcarriers
% CS : nominal constellation size
% Rb : bit rate (b/s)
% prefix : OFDM type {'SSB', 'DC', 'ACO'}
% power_allocation_type : {'Levin-Campello-MA', 'preemphasis'}
disp('-- SSB-OFDM simulation')
ofdm = ofdm(256, 208, 16, sim.Rb, 'SSB', 'Levin-Campello-MA');
ENOB = 5;

ofdm.set_cyclic_prefix(5, 5); % set cyclic prefix length. Should be consistent with channel memory length

%% Time and frequency
sim.N = sim.Mct*(ofdm.Nc + ofdm.Ncp)*sim.Nsymb;           % total number of points simulated in continuous time
sim.fs = ofdm.fs*sim.Mct;
[sim.f, sim.t] = freq_time(sim.N, sim.fs);

%% ==========================  Transmitter ================================
%% DAC
Tx.DAC.fs = ofdm.fs; % DAC sampling rate
Tx.DAC.ros = 1; % oversampling ratio of transmitter DSP
Tx.DAC.resolution = ENOB; % DAC effective resolution in bits
Tx.DAC.filt = design_filter('butter', 5, 0.5*ofdm.fs/(sim.fs/2)); % DAC analog frequency response

%% Modulator
Tx.rexdB = -15;  % extinction ratio in dB. Defined as Pmin/Pmax
Tx.rclip = 3; % clipping ratio
Tx.RIN = -150; % dB/Hz
Tx.alpha = 0;
Tx.CSPRdB = 15;

Tx.Mod.type = 'MZM';
Tx.Mod.Vswing = 0.6; % Voltage swing in MZM modulator normalized by Vpi/2
Tx.Mod.Vbias = 0.5; % Bias voltage normalized by Vpi/2
Tx.Mod.BW = 30e9; % moduator bandwidth
Tx.Mod.filt = design_filter('two-pole', Tx.Mod.BW, sim.fs);
% Tx.Mod.filt = design_filter('butter', 5, Tx.Mod.BW/(sim.fs/2));
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
% Constructor: OpticalAmplifier(Operation, param, Fn, Wavelength)
% - Opertation: either 'ConstantOutputPower' or 'ConstantGain'
% - param: GaindB if Operation = 'ConstantGain', or outputPower
% if Operation = 'ConstantOutputPower'
% - Fn:  noise figure in dB
% - Wavelength: operationl wavelength in m
% Rx.OptAmp = OpticalAmplifier('ConstantOutputPower', 0, 5, Tx.Laser.wavelength);
Rx.OptAmp = OpticalAmplifier('ConstantGain', 20, 5, Tx.Laser.wavelength);
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
%% ADC
Rx.ADC.ros = 1;
Rx.ADC.fs = ofdm.fs;
Rx.ADC.filt = design_filter('butter', 5, 0.5*ofdm.fs/(sim.fs/2)); % Antialiasing filter
Rx.ADC.ENOB = ENOB; % effective number of bits. Quantization is only included if sim.quantiz = true and ENOB ~= Inf
% Rx.ADC.rclip = 0.05;

%% Equalizer
Rx.AdEq.mu = 2e-3;
Rx.AdEq.Ntrain = 1224; % Number of frames used in training (if Inf all symbols are used)
% Note: training symbols are not used to compute the BER. Hence sim.Nsymb -
% Rx.AdEq.Ntrain must be large to obtain accurate BER estimate

%% Generate summary
ofdm_simulation_summary(sim, ofdm, Tx, Fibers, Rx);

%% Run simulation
[berOFDM, ofdm, OSNRdB] = ber_ssb_ofdm(ofdm, Tx, Fibers, Rx, sim)