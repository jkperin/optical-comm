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

%% Simulation parameters
sim.Nsymb = 2^15; % Number of symbols in montecarlo simulation
sim.Mct = 9;     % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
sim.L = 2;        % de Bruijin sub-sequence length (ISI symbol length)
sim.BERtarget = 1.8e-4; 
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

%
sim.OptimizeGain = true;
sim.shot = true; % include shot noise in montecarlo simulation (always included for pin and apd case)
sim.RIN = ~true; % include RIN noise in montecarlo simulation
sim.verbose = false; % show stuff

%% M-PAM
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
switch mpam.M
    case 4
        tx.PtxdBm = -30:1:-12;
    case 8
        tx.PtxdBm = -22:2:-4;
    case 16
       tx.PtxdBm = -18:2:-2;
end
   
tx.lamb = 1310e-9; % wavelength
tx.alpha = 0; % chirp parameter
tx.RIN = -145;  % dB/Hz
tx.rexdB = -10;  % extinction ratio in dB. Defined as Pmin/Pmax

% Modulator frequency response
tx.modulator.fc = 30e9; % modulator cut off frequency
tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
tx.modulator.h = @(t) [0*t(t < 0) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0))];
tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

%% Fiber
fiber = fiber(); % fiber(L, att(lamb), D(lamb))

%% Receiver
rx.N0 = (30e-12).^2; % thermal noise psd
rx.Id = 10e-9; % dark current
rx.R = 1; % responsivity
% Electric Lowpass Filter
% rx.elefilt = design_filter('bessel', 5, mpam.Rs/(sim.fs/2));
rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);
% rx.elefilt = design_filter('matched', @(t) conv(mpam.pshape(t), 1/sim.fs*tx.modulator.h(t/sim.fs), 'full') , 1/sim.Mct);

%% Equalization
rx.eq.type = 'Fixed TD-SR-LE';
% rx.eq.ros = 2;
rx.eq.Ntaps = 15;
% rx.eq.Ntrain = 2e3;
% rx.eq.mu = 1e-2;

%% PIN
% (GaindB, ka, GainBW, R, Id) 
pin = apd(0, 0, Inf, rx.R, rx.Id);

apd = apd(10, 0.1, 20e9, rx.R, rx.Id);

% Calculate BER
disp('BER with PIN')
sim.OptimizeGain = false;
ber_pin = apd_ber(mpam, tx, fiber, pin, rx, sim);

disp('BER with APD')
sim.OptimizeGain = true;
[ber_apd, ~, apd] = apd_ber(mpam, tx, fiber, apd, rx, sim);

varShot = apd.varShot(1e-3*10.^(tx.PtxdBm/10), rx.elefilt.noisebw(sim.fs/2));
varTherm = rx.N0*rx.elefilt.noisebw(sim.fs/2);

figure, box on, hold on
plot(tx.PtxdBm, 10*log10(varShot/1e-3), tx.PtxdBm([1 end]), 10*log10(varTherm/1e-3)*[1 1]);


%% Figures
figure, hold on, grid on, box on
plot(tx.PtxdBm, log10(ber_pin.gauss), '-k')
plot(tx.PtxdBm, log10(ber_apd.gauss), '-r')

plot(tx.PtxdBm, log10(ber_apd.count), '--or')
plot(tx.PtxdBm, log10(ber_pin.count), '--ok')

plot(tx.PtxdBm, log10(ber_apd.awgn), ':r')
plot(tx.PtxdBm, log10(ber_pin.awgn), ':k')

xlabel('Received Power (dBm)')
ylabel('log(BER)')
legend('PIN', 'APD', 'Location', 'SouthWest')
axis([tx.PtxdBm(1) tx.PtxdBm(end) -8 0])
set(gca, 'xtick', tx.PtxdBm)
    