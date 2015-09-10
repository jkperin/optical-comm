%% Power penalty vs APD bandwidth
clear, clc, close all


addpath ../f % general functions
addpath ../mpam/
addpath ../apd
addpath ../apd/f

%% Simulation parameters
sim.Nsymb = 2^15; % Number of symbols in montecarlo simulation
sim.Mct = 15;     % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
sim.L = 2;        % de Bruijin sub-sequence length (ISI symbol length)
sim.Me = 16; % number of used eigenvalues
sim.BERtarget = 1e-4; 
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

sim.polarizer = false;
sim.shot = true; % include shot noise in montecarlo simulation (always included for pin and apd case)
sim.RIN = true; % include RIN noise in montecarlo simulation
sim.verbose = false; % show stuff

%% M-PAM
mpam = PAM(4, 100e9, 'equally-spaced', @(n) double(n >= 0 & n < sim.Mct));

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
        tx.PtxdBm = -26:-8;
    case 8
        tx.PtxdBm = -22:2:-4;
    case 16
       tx.PtxdBm = -18:2:-2;
end
   
tx.lamb = 1310e-9; % wavelength
tx.alpha = 0; % chirp parameter
tx.RIN = -150;  % dB/Hz
tx.rexdB = -10;  % extinction ratio in dB. Defined as Pmin/Pmax

% Modulator frequency response
tx.modulator.fc = 20e9; % modulator cut off frequency
tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
tx.modulator.h = @(t) [0*t(t < 0) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0))];
tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

%% Fiber
fiber = fiber(); % fiber(L, att(lamb), D(lamb))
% fiber = fiber();

%% Receiver
rx.N0 = (30e-12).^2; % thermal noise psd
rx.Id = 10e-9; % dark current
rx.R = 1; % responsivity
% Electric Lowpass Filter
rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);

%% Equalization
rx.eq.type = 'Fixed TD-SR-LE';
% rx.eq.ros = 2;
rx.eq.Ntaps = 31;
% rx.eq.Ntrain = 2e3;
% rx.eq.mu = 1e-2;

%% PIN
% (GaindB, ka, BW, R, Id) 
pin = apd(0, 0, Inf, rx.R, rx.Id);


apd = apd(10, 0.1, 10e9, rx.R, rx.Id);

% 
Fc = (10:2.5:40)*1e9;

ber_pin = apd_ber(mpam, tx, fiber, pin, rx, sim);

Ptx_req_pin = interp1(ber_pin.gauss, tx.PtxdBm, sim.BERtarget);

% BER calculation
figure(1), hold on, grid on, box on
for k = 1:length(Fc)
    apd.BW = Fc(k);
    
    apd.Gain = apd.optGain(mpam, tx, fiber, rx, sim, 'margin');
    
    optGain(k) = apd.Gain;
    
    ber_apd(k) = apd_ber(mpam, tx, fiber, apd, rx, sim);
        
    Ptx_req_apd(k) = interp1(ber_apd(k).gauss, tx.PtxdBm, sim.BERtarget);
    
    %% Figures
    figure(1)
    plot(tx.PtxdBm, log10(ber_pin.gauss), '-b')
    plot(tx.PtxdBm, log10(ber_apd(k).gauss), '-r')
    plot(tx.PtxdBm, log10(ber_pin.awgn), '-k')
    plot(tx.PtxdBm, log10(ber_pin.awgn_ne), '--b')
    plot(tx.PtxdBm, log10(ber_apd(k).awgn_ne), '--r')

    plot(tx.PtxdBm, log10(ber_pin.count), '-ob')
    plot(tx.PtxdBm, log10(ber_apd(k).count), '-or')

    xlabel('Received Power (dBm)')
    ylabel('log(BER)')
    legend('Inf BW', 'Finite BW', 'Location', 'SouthWest')
    axis([tx.PtxdBm(1) tx.PtxdBm(end) -8 0])
    set(gca, 'xtick', tx.PtxdBm)
end

figure, hold on, box on
plot(Fc/1e9, Ptx_req_pin-Ptx_req_apd, '-o')
xlabel('Frequency (GHz)')
ylabel('Power margin improvement (dB)')


figure, hold on, box on
plot(Fc/1e9, optGain, '-o')
xlabel('Frequency (GHz)')
ylabel('Optimal gain')