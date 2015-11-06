%% Calculate BER of amplified IM-DD system. The only frequency responses taken 
% into account are the optical bandpass filter and the antialiasing electrical
% filter

% BER is calculated using the following methods
% 1. Montecarlo simulation
% 2. KLSE in frequency domain & Saddlepoint approximation of tail
% probabilities
% 3. KLSE with Fourier series & Saddlepoint approximation of tail
% probabilities
% 4. Gaussian approximation

clear, clc, close all

addpath ../mpam
addpath ../f % general functions
addpath f

% Simulation parameters
sim.Nsymb = 2^15; % Number of symbols in montecarlo simulation
sim.Mct = 15;      % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
sim.Me = 16;       % Number of used eigenvalues
sim.L = 2; % de Bruijin sub-sequence length (ISI symbol length)
sim.BERtarget = 1e-4; 
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

sim.polarizer = ~true;
sim.shot = true; % include shot noise in montecarlo simulation 
sim.RIN = true; % include RIN noise in montecarlo simulation
sim.verbose = false; % show stuff

% M-PAM
% M, Rb, leve_spacing, pshape
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
tx.PtxdBm = -24:1:-14;

tx.lamb = 1310e-9;
tx.RIN = -150;  % dB/Hz
tx.rexdB = -10;  % extinction ratio in dB. Defined as Pmin/Pmax

%% fiber
fiber = fiber(); % back-to-back

%% Receiver
rx.N0 = (30e-12).^2; % thermal noise psd
rx.Id = 10e-9; % dark current
rx.R = 1; % responsivity
% Electric Lowpass Filter
% rx.elefilt = design_filter('bessel', 5, 0.5*mpam.Rs/(sim.fs/2));
rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);
% Optical Bandpass Filter
rx.optfilt = design_filter('fbg', 4, 200e9/(sim.fs/2));

%% SOA
soa = soa(5, 7, 1310e-9, 20); % soa(GaindB, NF, lambda, maxGaindB)

%% Time and frequency
dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

sim.t = t;
sim.f = f;

% Calculate KLSE in the frequency domain
% [D_freq, Phi_freq, Fmax_freq] = klse_freq(rx, sim);

% KLSE Fourier Series Expansion 
[rx.U_fourier, rx.D_fourier, rx.Fmax_fourier] = klse_fourier(rx, sim, sim.Mct*(mpam.M^sim.L + 2*sim.L)); 

if sim.verbose
    figure, hold on
    plot(sim.f/1e9, abs(rx.elefilt.H(sim.f/sim.fs)).^2)
    plot(sim.f/1e9, abs(rx.optfilt.H(sim.f/sim.fs)).^2)
    legend('electrical', 'optical')
    xlabel('Frequency (GHz)')
    ylabel('|H(f)|^2')
    grid on
end

ber = soa_ber(mpam, tx, fiber, soa, rx, sim);

% ber_klse_freq = zeros(size(Ptx));
% for k = 1:length(Ptx)
%     tx.Ptx = Ptx(k);
% 
%     [ber_klse_freq(k), ber_pdf(k)] = ber_soa_klse_freq(D_freq, Phi_freq, Fmax_freq, mpam, tx, soa, rx, sim);
% end

figure(1), hold on, grid on
plot(tx.PtxdBm, log10(ber.count), '-o')
plot(tx.PtxdBm, log10(ber.est)) % KLSE Fourier
% plot(tx.PtxdBm, log10(ber_klse_freq))
plot(tx.PtxdBm, log10(ber.gauss))
plot(tx.PtxdBm, log10(ber.awgn))
legend('Monte Carlo', 'KLSE Fourier & Saddlepoint Approx',... %  'KLSE Frequency Domain & Saddlepoint Approx',...
        'Gaussian Approximation', 'AWGN approximation',...
    'Location', 'SouthWest')
axis([tx.PtxdBm([1 end]) -8 0])
xlabel('Transmitted Power (dBm)')
ylabel('log_{10}(BER)')
    