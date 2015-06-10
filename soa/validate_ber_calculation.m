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

addpath f

% Simulation parameters
sim.Nsymb = 2^16; % Number of symbols in montecarlo simulation
sim.Mct = 16;      % Oversampling ratio to simulate continuous time (must be even)  
sim.Me = 64; % Number of eigenvalues used only in klse_freq
sim.L = 2; % de Bruijin sub-sequence length (ISI symbol length)
sim.M = 2; % Ratio of optical filter BW and electric filter BW (must be integer)
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
PtxdBm = -25:-15;
Ptx = 1e-3*10.^(PtxdBm/10);
tx.rex = 10;  % extinction ratio in dB

%% Receiver
rx.N0 = (20e-12).^2; % thermal noise psd
% Electric Lowpass Filter
rx.elefilt = design_filter('bessel', 5, 1.25*mpam.Rs/(sim.fs/2));
% Optical Bandpass Filter
rx.optfilt = design_filter('fbg', 5, sim.M*rx.elefilt.fcnorm);

%% SOA
soa = soa(5, 9, 1310e-9, 20); % soa(GaindB, NF, lambda, maxGaindB)

%% Time and frequency
dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';
td = t(sim.Mct/2:sim.Mct:end);

sim.t = t;
sim.f = f;

% Generate unipolar PAM signal
if strcmp(mpam.level_spacing, 'uniform')
    mpam.a = (0:2:2*(mpam.M-1)).';
    mpam.b = (1:2:(2*(mpam.M-1)-1)).';
elseif strcmp(mpam.level_spacing, 'nonuniform')
   [mpam.a, mpam.b] = level_spacing_optm(mpam, tx, soa, rx, sim);
else
    error('Invalide Option!')
end

% Calculate KLSE in the frequency domain
[D_freq, Phi_freq, Fmax_freq] = klse_freq(rx, sim);

% KLSE Fourier Series Expansion 
[U_fourier, D_fourier, Fmax_fourier] = klse_fourier(rx, sim, sim.Mct*(mpam.M^sim.L + 2*sim.L)); 

if sim.verbose
    figure, hold on
    plot(sim.f/1e9, abs(rx.elefilt.H(sim.f/sim.fs)).^2)
    plot(sim.f/1e9, abs(rx.optfilt.H(sim.f/sim.fs)).^2)
    legend('electrical', 'optical')
    xlabel('Frequency (GHz)')
    ylabel('|H(f)|^2')
    grid on
end

ber_mc = zeros(size(Ptx));
ber_klse_fourier = zeros(size(Ptx));
ber_klse_freq = zeros(size(Ptx));
ber_gauss = zeros(size(Ptx));
for k = 1:length(Ptx)
    tx.Ptx = Ptx(k);

%     tx.Ptx = mean(mpam.a)/soa.Gain; % for testing level spacing optimization
    
    % Montecarlo simulation
    ber_mc(k) = ber_soa_montecarlo(mpam, tx, soa, rx, sim);
    
    % approx ber
%     [ber_klse_freq(k), ber_pdf(k)] = ber_soa_klse_freq(D_freq, Phi_freq, Fmax_freq, mpam, tx, soa, rx, sim);
%     [ber_klse_fourier(k), ber_gauss(k), ber_pdf(k)] = ber_soa_klse_fourier(U_fourier, D_fourier, Fmax_fourier, mpam, tx, soa, rx, sim);

    ber_klse_freq(k) = ber_soa_klse_freq(D_freq, Phi_freq, Fmax_freq, mpam, tx, soa, rx, sim);
    [ber_klse_fourier(k), ber_gauss(k)] = ber_soa_klse_fourier(U_fourier, D_fourier, Fmax_fourier, mpam, tx, soa, rx, sim);
    1;
end

figure, hold on
plot(PtxdBm, log10(ber_mc), '-o')
plot(PtxdBm, log10(ber_klse_fourier))
plot(PtxdBm, log10(ber_klse_freq))
plot(PtxdBm, log10(ber_gauss))
% plot(PtxdBm, log10(ber_pdf))
legend('Monte Carlo', 'KLSE Frequency Domain & Saddlepoint Approx', ...
    'KLSE Fourier & Saddlepoint Approx', 'Gaussian Approximation',...
    'Location', 'SouthWest')
axis([PtxdBm(1) PtxdBm(end) -10 0])
    