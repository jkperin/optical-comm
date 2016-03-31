%% Validate APD BER
clear, clc

addpath ../mpam
addpath ../f
addpath f

rx.eq.ros = 5/4; % oversampling ratio for DSP

% Simulation parameters
sim.Nsymb = 2^15; % Number of symbols in montecarlo simulation
sim.Mct = 8*rx.eq.ros;    % Oversampling ratio to simulate continuous time (must be odd so that sampling is done right, and FIR filters have interger grpdelay)  
sim.L = 4;        % de Bruijin sub-sequence length (ISI symbol length)
sim.BERtarget = 1.8e-4; 
sim.Ndiscard = 64;  % number of 0 symbols to be inserted at the begining and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation
sim.WhiteningFilter = true;

%
sim.RIN = true; % include RIN noise. Only included in montecarlo simulation

sim.plots = containers.Map();
sim.plots('BER') = 1;
sim.plots('Adaptation MSE') = 1;
sim.plots('Empirical noise pdf') = 0;
sim.plots('Frequency Response') = 0;

% M-PAM
mpam = PAM(4, 107e9, 'equally-spaced', @(n) double(n >= 0 & n < sim.Mct));

%% Time and frequency
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'

dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

sim.t = t;
sim.f = f;

%% Transmitter
tx.PtxdBm = -25:-15;

tx.lamb = 1250e-9; % wavelength
tx.alpha = 2; % chirp parameter
tx.RIN = -150;  % dB/Hz
tx.rexdB = -10;  % extinction ratio in dB. Defined as Pmin/Pmax

% Modulator frequency response
tx.modulator.BW = 30e9;
tx.modulator.fc = tx.modulator.BW/sqrt(sqrt(2)-1); % converts to relaxation frequency
tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
tx.modulator.h = @(t) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0));
tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

%% Fiber
fiber = fiber(1e3);

%% Receiver
rx.N0 = (30e-12).^2; % thermal noise psd
% Electric Lowpass Filter

%% Equalization
% rx.eq.type = 'None';
rx.eq.type = 'Adaptive TD-LE';
% rx.eq.type = 'Fixed TD-SR-LE';
rx.eq.Ntaps = 7;
rx.eq.mu = 1e-2;
rx.eq.Ntrain = Inf; % Number of symbols used in training (if Inf all symbols are used)
rx.eq.Ndiscard = [1e4 512]; % symbols to be discard from begining and end of sequence due to adaptation, filter length, etc
rx.elefilt = design_filter('butter', 5, rx.eq.ros/2*mpam.Rs/(sim.fs/2));
% rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);

%% APD 
% (GaindB, ka, [BW GBP=Inf], R, Id) 
apdG = apd(10, 0.2, [20e9 300e9], 1, 10e-9);
% apdG = apd(18, 0.3, Inf, 1, 10e-9);

% BER
sim.OptimizeGain = ~true;
ber_apd = apd_ber(mpam, tx, fiber, apdG, rx, sim);

mpam.level_spacing = 'optimized';
ber_apd_eq = apd_ber(mpam, tx, fiber, apdG, rx, sim);
        