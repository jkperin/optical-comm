%% Compare performance of M-PAM in 3 different scenarios
% 1. PIN receiver, no amplifier
% 2. SOA & PIN receiver
% 3. APD
clear, clc, close all

addpath ../f % general functions
addpath ../soa
addpath ../soa/f
addpath ../apd
addpath ../apd/f

% Simulation parameters
sim.Nsymb = 2^14; % Number of symbols in montecarlo simulation
sim.Mct = 8;     % Oversampling ratio to simulate continuous time (must be even)  
sim.L = 2;        % de Bruijin sub-sequence length (ISI symbol length)
sim.M = 4; % Ratio of optical filter BW and electric filter BW (must be integer)
sim.verbose = ~true; % show stuff
sim.BERtarget = 1e-4; 
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

sim.shot = true; % include shot noise in montecarlo simulation (always included for pin and apd case)
sim.RIN = true; % include RIN noise in montecarlo simulation

% M-PAM
mpam.level_spacing = 'uniform'; % M-PAM level spacing: 'uniform' or 'non-uniform'
mpam.M = 4;
mpam.Rb = 100e9;
mpam.Rs = mpam.Rb/log2(mpam.M);
mpam.pshape = @(n) ones(size(n)); % pulse shape

% 
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'

%% Transmitter
tx.PtxdBm = -30:2:-10;
tx.RIN = -150;  %dB/Hz
tx.rex = 10;  % extinction ratio in dB

%% Receiver
rx.N0 = (20e-12).^2; % thermal noise psd
rx.Id = 10e-9; % dark current
rx.R = 1; % responsivity
% Electric Lowpass Filter
% rx.elefilt = design_filter('bessel', 5, 1/2*mpam.Rs/(sim.fs/2));
rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);
% Optical Bandpass Filter
rx.optfilt = design_filter('fbg', 0, sim.M*rx.elefilt.fcnorm);

%% PIN
% (GaindB, ka, GainBW, R, Id) 
pin = apd(0, 0, Inf, rx.R, rx.Id);

%% APD 
% (GaindB, ka, GainBW, R, Id) 
% finite Gain x BW
apd_fin = apd(8.2761, 0.09, 340e9, rx.R, rx.Id); % gain optimized for uniformly-spaced 4-PAM with matched filter
% apd_fin.optimize_gain(mpam, tx, rx, sim);

if strcmp(mpam.level_spacing, 'uniform')
     % uniform, infinite Gain x BW (4-PAM)
    apd_inf = apd(11.0075, 0.09, Inf, 1, 10e-9); % gain optimized for 4-PAM with matched filter
%     apd_inf.optimize_gain(mpam, tx, rx, sim);

%     apd_inf = apd(8.4888, 0.09, Inf, rx.R, rx.Id); % uniform, infinite Gain x BW (8-PAM)
elseif strcmp(mpam.level_spacing, 'nonuniform')
    % nonuniform, infinite gain x BW
    apd_inf = apd(13.8408, 0.09, Inf, rx.R, rx.Id); % gain optimized for 4-PAM with matched filter
%     apd_inf.optimize_gain(mpam, tx, rx, sim);
end

%% SOA
% soa(GaindB, NF, lambda, maxGaindB)
soa = soa(20, 9, 1310e-9, 20); 

% BER
ber_soa = soa_ber(mpam, tx, soa, rx, sim);
ber_apd_fin = apd_ber(mpam, tx, apd_fin, rx, sim);
ber_apd_inf = apd_ber(mpam, tx, apd_inf, rx, sim);
ber_pin = apd_ber(mpam, tx, pin, rx, sim);

%% Figures
figure, hold on, grid on, box on
plot(tx.PtxdBm, log10(ber_soa.est), '-b')
plot(tx.PtxdBm, log10(ber_apd_fin.est), '--r')
plot(tx.PtxdBm, log10(ber_apd_inf.est), '--m')
plot(tx.PtxdBm, log10(ber_pin.gauss), '-k')

plot(tx.PtxdBm, log10(ber_soa.count), '-ob')
plot(tx.PtxdBm, log10(ber_apd_fin.count), '-or')
plot(tx.PtxdBm, log10(ber_apd_inf.count), '-om')
plot(tx.PtxdBm, log10(ber_pin.count), '-ok')

plot(tx.PtxdBm, log10(ber_soa.gauss), '--b')
plot(tx.PtxdBm, log10(ber_apd_fin.gauss), '-r')
plot(tx.PtxdBm, log10(ber_apd_inf.gauss), '-m')

xlabel('Received Power (dBm)')
ylabel('log(BER)')
legend('SOA', 'APD Gain x BW = 340 GHz', 'APD Gain x BW = Inf', 'PIN', 'Location', 'SouthWest')
axis([tx.PtxdBm(1) tx.PtxdBm(end) -8 0])
set(gca, 'xtick', tx.PtxdBm)
    