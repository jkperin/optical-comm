%% Validate APD BER
clear, clc

addpath ../mpam
addpath ../f
addpath f

rx.eq.ros = 5/4;

% Simulation parameters
sim.Nsymb = 2^18; % Number of symbols in montecarlo simulation
sim.Mct = 8*rx.eq.ros;    % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
sim.L = 4;        % de Bruijin sub-sequence length (ISI symbol length)
sim.BERtarget = 1.8e-4; 
sim.Ndiscard = 64; % number of symbols to be discarded from the begining and end of the sequence (should be larger than Ntrain, Ntaps, Nfft, etc
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation
sim.WhiteningFilter = true;

%
sim.shot = true; % include shot noise. Only included in montecarlo simulation (except for APD)
sim.RIN = true; % include RIN noise. Only included in montecarlo simulation
sim.verbose = 1; % verbose level: verbose is decremented on each function. If verbose=0, nothing is shown

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
tx.PtxdBm = -25:-10;

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
fiber = fiber();

%% Receiver
rx.N0 = (30e-12).^2; % thermal noise psd
% Electric Lowpass Filter

%% Equalization
% rx.eq.type = 'None';
rx.eq.type = 'Adaptive TD-LE';
rx.eq.Ntaps = 11;
rx.eq.mu = 1e-2;
rx.eq.Ntrain = Inf; % Number of symbols used in training (if Inf all symbols are used)
rx.eq.Ndiscard = [4e4 64]; % symbols to be discard from begining and end of sequence due to adaptation, filter length, etc
rx.elefilt = design_filter('butter', 5, rx.eq.ros/2*mpam.Rs/(sim.fs/2));
% rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);

% rx.eq.type = 'Fixed FD-FS-LE';
% rx.eq.ros = 2;
% rx.eq.Ntaps = 31;
% rx.eq.Nfft = 512;
% rx.eq.Noverlap = 256;
% rx.eq.Ntrain = 2e3;
% rx.eq.mu = 1e-2;

%% APD 
% (GaindB, ka, [BW GBP=Inf], R, Id) 
apdG = apd(12, 0.1, [20e9 300e9], 1, 10e-9);

% BER
sim.OptimizeGain = ~true;
ber_apd = apd_ber(mpam, tx, fiber, apdG, rx, sim);

% mpam.level_spacing = 'optimized';
% ber_apd_eq = apd_ber(mpam, tx, fiber, apdG, rx, sim);
        