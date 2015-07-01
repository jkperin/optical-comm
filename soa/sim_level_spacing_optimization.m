%% Simulations for SOA with nonuniform and uniform level spacing compared with Gaussian approximation
clear, clc, close all

addpath ../f % general functions
addpath ../soa
addpath ../soa/f
addpath ../apd
addpath ../apd/f

% Simulation parameters
sim.Nsymb = 2^15; % Number of symbols in montecarlo simulation
sim.Mct = 16;     % Oversampling ratio to simulate continuous time (must be even)  
sim.L = 2;        % de Bruijin sub-sequence length (ISI symbol length)
sim.M = 4; % Ratio of optical filter BW and electric filter BW (must be integer)
sim.Me = 16; % number of used eigenvalues
sim.verbose = ~true; % show stuff
sim.BERtarget = 1e-4; 
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

sim.shot = false; % include shot noise in montecarlo simulation (always included for pin and apd case)
sim.RIN = false; % include RIN noise in montecarlo simulation

% M-PAM
mpam.level_spacing = 'uniform'; % M-PAM level spacing: 'uniform' or 'non-uniform'
mpam.M = 4;
mpam.Rb = 100e9;
mpam.Rs = mpam.Rb/log2(mpam.M);
mpam.pshape = @(n) double(n >= 0 & n < sim.Mct); % pulse shape

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
tx.alpha = 0; % chirp parameter
tx.RIN = -150;  % dB/Hz
tx.rexdB = -15;  % extinction ratio in dB. Defined as Pmin/Pmax

% Modulator frequency response
tx.kappa = 1; % controls attenuation of I to P convertion

%% Fiber
fiber = fiber(); % fiber(L, att(lamb), D(lamb))

%% Receiver
rx.N0 = (20e-12).^2; % thermal noise psd
rx.Id = 10e-9; % dark current
rx.R = 1; % responsivity
% Electric Lowpass Filter
rx.elefilt = design_filter('bessel', 5, mpam.Rs/(sim.fs/2));
% rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);
% Optical Bandpass Filter
rx.optfilt = design_filter('butter', 4, sim.M*rx.elefilt.fcnorm);

% KLSE Fourier Series Expansion (done here because depends only on filters
% frequency response)
% klse_fourier(rx, sim, N, Hdisp)
[rx.U_fourier, rx.D_fourier, rx.Fmax_fourier] = klse_fourier(rx, sim, sim.Mct*(mpam.M^sim.L + 2*sim.L)); 

%% SOA
% soa(GaindB, NF, lambda, maxGaindB)
soa = soa(10, 9, 1310e-9, 20); 

FndB = 3:10;

for k = 1:length(FndB)
    soa.Fn = FndB(k);
    
    mpam.level_spacing = 'uniform'; % M-PAM level spacing: 'uniform' or 'non-uniform'
    
    [ber_uniform(k)] = soa_ber(mpam, tx, fiber, soa, rx, sim);
    
    Psensitivity_uniform(k) = interp1(log10(ber_uniform(k).est), tx.PtxdBm, log10(sim.BERtarget), 'spline');
    
    mpam.level_spacing = 'nonuniform'; % M-PAM level spacing: 'uniform' or 'non-uniform'
    
    [ber_nonuniform(k)] = soa_ber(mpam, tx, fiber, soa, rx, sim);
    
    Psensitivity_nonuniform(k) = interp1(log10(ber_nonuniform(k).est), tx.PtxdBm, log10(sim.BERtarget), 'spline');  
    
    figure(1), hold on, grid on, box on
    plot(tx.PtxdBm, log10(ber_uniform(k).est), '-b')
    plot(tx.PtxdBm, log10(ber_uniform(k).count), ':ob')
    plot(tx.PtxdBm, log10(ber_uniform(k).gauss), '--b')  
    
    plot(tx.PtxdBm, log10(ber_nonuniform(k).est), '-r')
    plot(tx.PtxdBm, log10(ber_nonuniform(k).count), ':or')
    plot(tx.PtxdBm, log10(ber_nonuniform(k).gauss), '--r')
end

xlabel('Received Power (dBm)')
ylabel('log(BER)')
legend('KLSE Fourier', 'Montecarlo', 'Gaussian Approximmation', 'Location', 'SouthWest')
axis([tx.PtxdBm(1) tx.PtxdBm(end) -8 0])
set(gca, 'xtick', tx.PtxdBm)

%% Figures
figure, hold on, grid on, box on
plot(FndB, Psensitivity_uniform)
plot(FndB, Psensitivity_nonuniform)
legend('Uniform Level Spacing', 'Non-Uniform Level Spacing', 'Location', 'NorthWest')
xlabel('Noise Figure (dB)')
ylabel('Receiver Sensitivity (dBm)')


%% Plot Frequency response
% signal = design_filter('matched', mpam.pshape, 1/sim.Mct);
% Hsig = signal.H(sim.f/sim.fs); % signal frequency response
% figure, box on, grid on, hold on
% plot(f/1e9, abs(Hsig).^2)
% if isfield(tx, 'modulator')
%     plot(f/1e9, abs(tx.modulator.H(f)).^2)
% else
%     plot(f/1e9, ones(size(f)))
% end
% plot(f/1e9, abs(fiber.Hfiber(f, tx)).^2)
% plot(f/1e9, abs(rx.optfilt.H(f/sim.fs)).^2)
% plot(f/1e9, abs(rx.elefilt.H(f/sim.fs)).^2)
% legend('Signal', 'Modulator', 'Fiber frequency response (small-signal)', 'Optical filter', 'Receiver electric filter')
% xlabel('Frequency (GHz)')
% ylabel('|H(f)|^2')
% axis([0 rx.Fmax_fourier*sim.fs/1e9 0 3])
    