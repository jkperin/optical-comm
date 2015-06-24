%% Validate APD Gain Optimization
%% !! Needs equalization to work with matched filter or other filters that distort the signal
clear, clc, close all

addpath ../f
addpath f

% Simulation parameters
sim.Nsymb = 2^18; % Number of symbols in montecarlo simulation
sim.Mct = 16;     % Oversampling ratio to simulate continuous time (must be even)  
sim.L = 3;        % de Bruijin sub-sequence length (ISI symbol length)
sim.shot = true; % include shot noise. Only included in montecarlo simulation (except for APD)
sim.rin = ~true; % include RIN noise. Only included in montecarlo simulation
sim.verbose = ~true; % show stuff
sim.BERtarget = 1e-6; 
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

% M-PAM
mpam.level_spacing = 'nonuniform'; % M-PAM level spacing: 'uniform' or 'non-uniform'
mpam.M = 4;
mpam.Rb = 100e9;
mpam.Rs = mpam.Rb/log2(mpam.M);
mpam.pshape = @(n) ones(size(n)); % pulse shape

%% Time and frequency
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'

dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

sim.t = t;
sim.f = f;

%% Transmitter
tx.lamb = 1310e-9; % wavelength
tx.alpha = 2; % chirp parameter
tx.RIN = -150;  % dB/Hz
tx.rexdB = -Inf;  % extinction ratio in dB. Defined as Pmin/Pmax

% Modulator frequency response
tx.kappa = 1; % controls attenuation of I to P convertion
% tx.modulator.fc = 2*mpam.Rs; % modulator cut off frequency
% tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
% tx.modulator.h = @(t) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0));
% tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

%% Fiber
fiber = fiber();

%% Receiver
rx.N0 = (20e-12).^2; % thermal noise psd
% Electric Lowpass Filter
% rx.elefilt = design_filter('bessel', 5, 1*mpam.Rs/(sim.fs/2));
rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);

%% APD 
% (GaindB, ka, GainBW, R, Id)  
apd = apd(10.0851, 0.09, Inf, 1, 10e-9); % uniform, infinite gain x BW product

% Level spacing and decision threshold at receiver to achieve sim.BERtarget
[mpam.a, mpam.b] = level_spacing_optm_gauss_approx(mpam, tx, apd, rx, sim);

% Refer to transmitter
link_gain = tx.kappa*fiber.link_attenuation(tx.lamb)*apd.R*apd.Gain;
mpam.a = mpam.a/link_gain;
mpam.b = mpam.b/link_gain;

tx.Ptx = mean(mpam.a);

ber_count = ber_apd_montecarlo(mpam, tx, fiber, apd, rx, sim);

[~, ber_gauss] = ber_apd_doubly_stochastic(mpam, tx, fiber, apd, rx, sim);

[sim.BERtarget, ber_count, ber_gauss]

    