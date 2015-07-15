%% Vary Gain of Optical Amplifier

clear, clc, close all

addpath ../f % general functions
addpath f

profile on

% Simulation parameters
sim.Nsymb = 2^14; % Number of symbols in montecarlo simulation
sim.Mct = 17;     % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
sim.L = 3;        % de Bruijin sub-sequence length (ISI symbol length)
sim.Me = 8; % Number of used eigenvalues
sim.BERtarget = 1e-4; 
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

sim.shot = false; % include shot noise in montecarlo simulation 
sim.RIN = false; % include RIN noise in montecarlo simulation
sim.verbose = false; % show stuff

% M-PAM
mpam.level_spacing = 'uniform'; % M-PAM level spacing: 'uniform' or 'non-uniform'
mpam.M = 8;
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
if mpam.M == 8
    tx.PtxdBm = -22:2:-4;
else
    tx.PtxdBm = -26:1:-10;
end

tx.lamb = 1310e-9; % wavelength
tx.alpha = 2; % chirp parameter
tx.RIN = -150;  % dB/Hz
tx.rexdB = -Inf;  % extinction ratio in dB. Defined as Pmin/Pmax

tx.kappa = 1;
% tx.modulator.fc = 2*mpam.Rs; % modulator cut off frequency
% tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
% tx.modulator.h = @(t) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0));
% tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

%% Receiver
rx.N0 = (20e-12).^2; % thermal noise psd
rx.Id = 10e-9; % dark current
rx.R = 1; % responsivity
% Electric Lowpass Filter
% rx.elefilt = design_filter('bessel', 5, 0.5*mpam.Rs/(sim.fs/2));
rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);
% Optical Bandpass Filter
rx.optfilt = design_filter('fbg', 0, 200e9/(sim.fs/2));

% KLSE Fourier Series Expansion (done here because depends only on filters
% frequency response)
[rx.U_fourier, rx.D_fourier, rx.Fmax_fourier] = klse_fourier(rx, sim, sim.Mct*(mpam.M^sim.L + 2*sim.L)); 
%% SOA
% soa(GaindB, NF, lambda, maxGaindB)
soa_opt = soa(10, 9, tx.lamb, Inf);
soa_opt.optimize_gain(mpam, tx, fiber, rx, sim)

profile off
profile viewer

GainsdB = [5:2.5:20 soa_opt.GaindB];

soaG = soa(10, 9, 1310e-9, 20); 

figure, hold on, grid on
legends = {};
for k= 1:length(GainsdB)

    % Selected gain
    soaG.GaindB = GainsdB(k);

    % BER
    ber_soa = soa_ber(mpam, tx, fiber, soaG, rx, sim);
    
%     plot(tx.PtxdBm, log10(ber_soa.count), '-o')
    plot(tx.PtxdBm, log10(ber_soa.est), '-')
%     plot(tx.PtxdBm, log10(ber_soa.gauss), '--')
    legends = [legends, sprintf('Gain = %.1f dB', GainsdB(k))];
end
xlabel('Received Power (dBm)')
ylabel('log(BER)')
legend(legends{:})





    