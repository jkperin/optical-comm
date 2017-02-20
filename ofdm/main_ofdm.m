%% Evaluation of OFDM in IM-DD system, which may be amplified or not
clear, clc, close all

addpath f/
addpath ../f/
addpath ../apd/

%% Transmit power swipe
% Tx.PtxdBm = -10:-5; % transmitter power range
Tx.PtxdBm = -12:-7; % transmitter power range

%% Simulation parameters
sim.Rb = 112e9;    % bit rate in bits/sec
sim.Nsymb = 2^11; % Number of symbols in montecarlo simulation
sim.Mct = 5;      % Oversampling ratio to simulate continuous time. Must be integer multiple of sim.ros.txDSP and numerator of sim.ros.rxDSP
sim.BERtarget = 1.8e-4; 
sim.Ndiscard = 128; % number of symbols to be discarded from the begining and end of the sequence
sim.Modulator = 'DML'; % 'MZM' or 'DML'
sim.OFDM = 'ACO-OFDM'; % {'DC-OFDM', 'ACO-OFDM'}
 
%% Simulation control
sim.preAmp = ~true;
sim.RIN = true; % include RIN noise in montecarlo simulation
sim.phase_noise = true; % whether to simulate laser phase noise
sim.PMD = false; % whether to simulate PMD
sim.quantiz = true; % whether quantization is included
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
sim.Plots('Decision errors') = 0;
sim.Plots('Channel frequency response') = 0;
sim.Plots('OSNR') = 0;
sim.shouldPlot = @(x) sim.Plots.isKey(x) && sim.Plots(x);

%% OFDM 
% OFDM constructor: ofdm(Nc, Nu, CS, Rb, power_allocation_type (optional, default = 'palloc')
% Nc : number of subcarriers
% Nu : number of used subcarriers
% CS : nominal constellation size
% Rb : bit rate (b/s)
% power_allocation_type : {'palloc', 'preemphasis'}
if strcmpi(sim.OFDM, 'ACO-OFDM')
    disp('-- ACO-OFDM simulation')
    ofdm = ofdm(256, 208, 64, sim.Rb, 'palloc'); 
    ofdm.aco_ofdm_config();
    ENOB = 6;
elseif strcmpi(sim.OFDM, 'DC-OFDM')
    disp('-- DC-OFDM simulation')
    ofdm = ofdm(256, 208, 16, sim.Rb, 'palloc');
    ENOB = 5;
else
    error('sim.OFDM must be either DC-OFDM or ACO-OFDM')
end
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
Tx.rclip = 4.5; % clipping ratio
Tx.alpha = 0; % chirp parameter
Tx.RIN = -150; % dB/Hz
Tx.Mod.type = sim.Modulator;   
Tx.Mod.Vswing = 1; % voltage swing in MZM modulator. Thi is normalized by Vpi. Vswing = 1 corresponds to swing from 0 to Vpi.
Tx.Mod.BW = 50e9; % moduator bandwidth
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
Rx.AdEq.mu = 1e-3;
Rx.AdEq.Ntrain = 1224; % Number of frames used in training (if Inf all symbols are used)
% Note: training symbols are not used to compute the BER. Hence sim.Nsymb -
% Rx.AdEq.Ntrain must be large to obtain accurate BER estimate

%% Generate summary
ofdm_simulation_summary(sim, ofdm, Tx, Fibers, Rx);

%% Run simulation
[berOFDM, ofdm, OSNRdB] = ber_ofdm(ofdm, Tx, Fibers, Rx, sim)