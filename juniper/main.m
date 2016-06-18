%% Calculate BER of amplified IM-DD system. 
% - CD is partially compensated by DCF
% - Equalization is done using a fractionally spaced linear equalizer
% Simulations include modulator (as a 2nd-order filter), fiber, optical 
% amplifier characterized only by gain and noise figure, optical bandpass 
% filter, antialiasing filter, sampling, and linear equalization

clear, clc

addpath f/ % Juniper project specific functions
addpath ../mpam % PAM
addpath ../f % general functions
addpath ../soa % for pre-amplifier 
addpath ../apd % for PIN photodetectors

rx.eq.ros = 5/4; % oversampling ratio for DSP
% if equalization type is fixed time-domain equalizer, then ros = 1

% Simulation parameters
sim.Rb = 56e9;    % bit rate in bits/sec
sim.Nsymb = 2^15; % Number of symbols in montecarlo simulation
sim.Mct = 15;      % Oversampling ratio to simulate continuous time
sim.Me = 16;       % Number of used eigenvalues
sim.L = 2; % de Bruijin sub-sequence length (ISI symbol length)
sim.BERtarget = 1e-4; 
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation
 
sim.RIN = true; % include RIN noise in montecarlo simulation
sim.stopSimWhenBERreaches0 = true; % stop simulation when counted BER reaches 0

% Plot
sim.Plots = containers.Map();
sim.Plots('BER') = 1;
sim.Plots('Equalizer') = 0;
sim.Plots('Adaptation MSE') = 0;
sim.Plots('Ouptut signal') = 0;
sim.Plots('Heuristic noise pdf') = 0;
sim.Plots('Frequency Response') = 0;

% M-PAM
% PAM(M, bit rate, leve spacing : {'equally-spaced', 'optimized'}, pulse
% shape: anonymous function containing pulse shape versus time samples)
mpam = PAM(4, sim.Rb, 'equally-spaced', @(n) double(n >= 0 & n < sim.Mct));

%% Time and frequency
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'

dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt);
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df);

sim.t = t;
sim.f = f;

%% Transmitter
tx.PtxdBm = -30:-15; % transmitter power range
% tx.PtxdBm = -15; % transmitter power range

tx.lamb = 1550e-9; % transmitter wavelength
tx.RIN = -150;  % dB/Hz
tx.rexdB = -10;  % extinction ratio in dB. Defined as Pmin/Pmax
tx.alpha = 0; % chirp parameter

% Modulator frequency response
tx.modulator.BW = 30e9;
tx.modulator.fc = tx.modulator.BW/sqrt(sqrt(2)-1); % converts to relaxation frequency
tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
tx.modulator.h = @(t) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0));
tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

%% Fiber
% fiber(Length in m, anonymous function for attenuation versus wavelength
% (default: att(lamb) = 0 i.e., no attenuation), anonymous function for 
% dispersion versus wavelength (default: SMF28 with lamb0 = 1310nm, S0 = 0.092 s/(nm^2.km))
SMF = fiber(80e3); % SMF fiber
DCF = fiber(32.5e3, @(lamb) 0, @(lamb) -0.1*(lamb-1550e-9)*1e3 - 40e-6); % DCF fiber

Fibers = [SMF DCF];

%% Amplifier
% Class SOA characterizes amplifier in terms of gain and noise figure
EDFA = soa(20, 5, 1550e-9, 20); % soa(GaindB, NFdB, lambda, maxGaindB)

%% Receiver
% pin(R: Responsivity (A/W), Id: Dark current (A), BW: bandwidth (Hz))
% PIN frequency response is modeled as a first-order system
rx.PD = pin(1, 10e-9, Inf);

% One-sided thermal noise PSD
rx.N0 = (30e-12).^2; 

% Optical Bandpass Filter (fiber Brag gratting)
rx.optfilt = design_filter('fbg', 4, 200e9/(sim.fs/2));

%% Equalizer
% TD = time domain, SR = symbol-rate, LE = linear equalizer

% rx.eq.type = 'None';
rx.eq.type = 'Adaptive TD-LE';
% rx.eq.type = 'Fixed TD-SR-LE';
rx.eq.Ntaps = 3;
rx.eq.mu = 1e-2;
rx.eq.Ntrain = 0.5e4; % Number of symbols used in training (if Inf all symbols are used)
rx.eq.Ndiscard = [1e4 512]; % symbols to be discard from begining and end of sequence due to adaptation, filter length, etc
rx.elefilt = design_filter('butter', 5, rx.eq.ros/2*mpam.Rs/(sim.fs/2));
% rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);

% % Calculate KLSE in the frequency domain
% % [D_freq, Phi_freq, Fmax_freq] = klse_freq(rx, sim);
% % 
% % KLSE Fourier Series Expansion 
% % [rx.U_fourier, rx.D_fourier, rx.Fmax_fourier] = klse_fourier(rx, sim, sim.Mct*(mpam.M^sim.L + 2*sim.L)); 
% % 
% % if sim.verbose
% %     figure, hold on
% %     plot(sim.f/1e9, abs(rx.elefilt.H(sim.f/sim.fs)).^2)
% %     plot(sim.f/1e9, abs(rx.optfilt.H(sim.f/sim.fs)).^2)
% %     legend('electrical', 'optical')
% %     xlabel('Frequency (GHz)')
% %     ylabel('|H(f)|^2')
% %     grid on
% % end

%% Generate summary
[simTable, txTable, rxTable] = generate_summary(sim, tx, rx);
SMF.summary(tx.lamb);
DCF.summary(tx.lamb);
EDFA.summary();
rx.PD.summary();

%% Run simulation
ber = preamplified_sys_ber(mpam, tx, Fibers, EDFA, rx, sim);