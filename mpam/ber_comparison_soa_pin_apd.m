%% Compare performance of M-PAM in 3 different scenarios
% 1. PIN receiver, no amplifier
% 2. SOA & PIN receiver
% 3. APD
clear, clc, close all

addpath f
addpath ../soa
addpath ../apd

% Simulation parameters
sim.Nsymb = 2^16; % Number of symbols in montecarlo simulation
sim.Mct = 16;     % Oversampling ratio to simulate continuous time (must be even)  
sim.L = 2;        % de Bruijin sub-sequence length (ISI symbol length)
sim.M = 2; % Ratio of optical filter BW and electric filter BW (must be integer)
sim.shot = true; % include shot noise. Only included in montecarlo simulation (except for APD)
sim.rin = true; % include RIN noise. Only included in montecarlo simulation
sim.verbose = ~true; % show stuff
sim.BERtarget = 1e-4; 
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

% M-PAM
mpam.level_spacing = 'uniform'; % M-PAM level spacing: 'uniform' or 'non-uniform'
mpam.M = 4;
mpam.Rb = 100e9;
mpam.Rs = mpam.Rb/log2(mpam.M);
mpam.pshape = @(n) ones(size(n)); % pulse shape

% 
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'

%% Transmitter
tx.PtxdBm = -25:-15;
tx.rex = 10;  % extinction ratio in dB

%% Receiver
rx.N0 = (20e-12).^2; % thermal noise psd
% Electric Lowpass Filter
rx.elefilt = design_filter('bessel', 5, 1.25*mpam.Rs/(sim.fs/2));
% Optical Bandpass Filter
rx.optfilt = design_filter('fbg', 5, sim.M*rx.elefilt.fcnorm);

%% PIN
% (GainBW, ka, Gain, noise_stats, R, Id) 
pin = apd(Inf, 0, 1, 'Gaussian', 1, 10e-9);

%% APD 
% (GainBW, ka, Gain, noise_stats, R, Id) 
apd = apd(Inf, 0.09, 10, 'Gaussian', 1, 10e-9);

%% SOA
% soa(GaindB, NF, lambda, maxGaindB)
soa = soa(5, 9, 1310e-9, 20); 

% BER for SOA system
ber_soa = soa_ber(mpam, tx, soa, rx, sim);

%% Figures
figure, hold on
plot(tx.PtxdBm, log10(ber_soa.count), '-o')
plot(tx.PtxdBm, log10(ber_soa.est))
plot(tx.PtxdBm, log10(ber_soa.gauss))
legend('Counted', 'KLSE Fourier & Saddlepoint Approx', 'Gaussian Approximation',...
    'Location', 'SouthWest')
axis([tx.PtxdBm(1) tx.PtxdBm(end) -10 0])
    