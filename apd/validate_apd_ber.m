%% Validate APD BER
clear, clc, close all

addpath ../mpam
addpath ../f
addpath f

% Simulation parameters
sim.Nsymb = 2^16; % Number of symbols in montecarlo simulation
sim.Mct = 15;    % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
sim.L = 2;        % de Bruijin sub-sequence length (ISI symbol length)
sim.BERtarget = 1.8e-4; 
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

%
sim.shot = true; % include shot noise. Only included in montecarlo simulation (except for APD)
sim.RIN = false; % include RIN noise. Only included in montecarlo simulation
sim.verbose = ~true; % show stuff

% M-PAM
mpam = PAM(4, 100e9, 'optimized', @(n) double(n >= 0 & n < sim.Mct));

%% Time and frequency
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'

dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

sim.t = t;
sim.f = f;

%% Transmitter
tx.PtxdBm = -20:-10;

tx.lamb = 1310e-9; % wavelength
tx.alpha = 0; % chirp parameter
tx.RIN = -150;  % dB/Hz
tx.rexdB = -10;  % extinction ratio in dB. Defined as Pmin/Pmax

% Modulator frequency response
tx.modulator.fc = 30e9; % modulator cut off frequency
tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
tx.modulator.h = @(t) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0));
tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

%% Fiber
fiber = fiber();

%% Receiver
rx.N0 = (30e-12).^2; % thermal noise psd
% Electric Lowpass Filter
% rx.elefilt = design_filter('bessel', 5, 1.25*mpam.Rs/(sim.fs/2));
% rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);

%% Equalization
% rx.eq.type = 'None';
rx.eq.type = 'Fixed TD-SR-LE';
% rx.eq.ros = 2;
rx.eq.Ntaps = 30;
% rx.eq.Ntrain = 2e3;
% rx.eq.mu = 1e-2;

%% APD 
% (GaindB, ka, BW, R, Id) 
apdG = apd(20, 0.1, Inf, 1, 10e-9);

% BER
sim.OptimizeGain = true;
ber_apd = apd_ber(mpam, tx, fiber, apdG, rx, sim);
     
figure, hold on, box on
plot(tx.PtxdBm, log10(ber_apd.gauss), '-');
plot(tx.PtxdBm, log10(ber_apd.count), '--o')
plot(tx.PtxdBm, log10(ber_apd.awgn), '--')    
legend('Estimated', 'Montecarlo', 'AWGN')
xlabel('Received Power (dBm)')
ylabel('log_{10}(BER)')
grid on
axis([tx.PtxdBm([1 end]) -8 0])


    