%% Calculate BER of amplified IM-DD system. The only frequency responses taken 
% into account are the optical bandpass filter and the antialiasing electrical
% filter
% BER is calculated via montecarlo simulation and through saddlepoint approximation
% to calculate tail probabilities given the moment generating funcion.
% Only ASE and thermal noise are considered
clear, clc, close all

addpath f

% Simulation parameters
sim.Nsymb = 2^18; % Number of symbols in montecarlo simulation
sim.Mct = 8;      % Oversampling ratio to simulate continuous time  
sim.Me = 32; % Number of eigenvalues of KL expansion
sim.L = 4; % de Bruijin sub-sequence length (ISI symbol length)
sim.M = 1; % Ratio of optical filter BW and electric filter BW (must be integer)
sim.levelSpacing = 'uniform'; % M-PAM level spacing ## non-uniform to be implemented
sim.verbose = ~true; % show stuff
sim.BERtarget = 1e-4; 
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

% M-PAM
mpam.M = 4;
mpam.Rb = 100e9;
mpam.Rs = mpam.Rb/log2(mpam.M);
mpam.pshape = @(n) ones(size(n)); % pulse shape

% 
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'

% Generate unipolar PAM signal
if strcmp(sim.levelSpacing, 'uniform')
    mpam.a = (0:2:2*(mpam.M-1)).';
    mpam.b = (1:2:(2*(mpam.M-1)-1)).';
elseif strcmp(sim.levelSpacing, 'nonuniform')
   %% ## to be implemented
else
    error('Invalide Option!')
end

PtxdBm = -20:-10;
Ptx = 1e-3*10.^(PtxdBm/10);
tx.rex = 10;  % extinction ratio in dB

%% Receiver
rx.N0 = (30e-12).^2;

% Electric Lowpass Filter
rx.elefilt = design_filter('bessel', 5, 1.25*mpam.Rs/(sim.fs/2));

% Optical Bandpass Filter
rx.optfilt = design_filter('butter', 5, sim.M*rx.elefilt.fcnorm);

if rx.optfilt.fcnorm > 0.5
    warning(sprintf('Oversampling ratio for continuous time simulation might not be high enough to model optical bandpass filter accurately: f3dB = %g GHz, fs = %g GHz', rx.optfilt.fcnorm*sim.fs/2e9, sim.fs/1e9));
end

% SOA
soa = soa(10^(10/10), 9, 1310e-9, 20);

% Time and frequency
dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';
td = t(sim.Mct/2:sim.Mct:end);

sim.t = t;
sim.f = f;

% Calculate KL series expansion in the frequency domain
[D, Phi, Fmax] = klse_freq(rx, sim);

if ~sim.verbose
    figure, hold on
    plot(sim.f/1e9, abs(rx.elefilt.H(sim.f/sim.fs)).^2)
    plot(sim.f/1e9, abs(rx.optfilt.H(sim.f/sim.fs)).^2)
    legend('electrical', 'optical')
    xlabel('Frequency (GHz)')
    ylabel('|H(f)|^2')
    grid on
end

ber_mc = zeros(size(Ptx));
ber_tail = zeros(size(Ptx));
for k = 1:length(Ptx)
    tx.Ptx = Ptx(k);
    
    % Montecarlo simulation
    ber_mc(k) = ber_soa_montecarlo(mpam, tx, soa, rx, sim);
    
    % approx ber
    ber_tail(k) = ber_soa_klse_freq(D, Phi, Fmax, mpam, tx, soa, rx, sim);
    1;
end

figure
plot(PtxdBm, log10(ber_mc), '-o', PtxdBm, log10(ber_tail))
legend('Monte Carlo', 'Saddlepoint Approx')
axis([PtxdBm(1) PtxdBm(end) -10 0])
    