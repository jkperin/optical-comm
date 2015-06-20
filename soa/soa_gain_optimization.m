%% Vary Gain of Optical Amplifier

clear, clc, close all

addpath ../f % general functions
addpath f

% Simulation parameters
sim.Nsymb = 2^14; % Number of symbols in montecarlo simulation
sim.Mct = 8;     % Oversampling ratio to simulate continuous time (must be even)  
sim.L = 2;        % de Bruijin sub-sequence length (ISI symbol length)
sim.M = 4; % Ratio of optical filter BW and electric filter BW (must be integer)
sim.Me = 16; % Number of used eigenvalues
sim.verbose = ~true; % show stuff
sim.BERtarget = 1e-4; 
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation
sim.shot = true; % include shot noise in montecarlo simulation 
sim.RIN = true; % include RIN noise in montecarlo simulation

% M-PAM
mpam.level_spacing = 'uniform'; % M-PAM level spacing: 'uniform' or 'non-uniform'
mpam.M = 4;
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
tx.PtxdBm = -25:-10;
tx.rex = 10;  % extinction ratio in dB
tx.RIN = -150;  % dB/Hz
tx.rex = 10;  % extinction ratio in dB

%% Receiver
rx.N0 = (20e-12).^2; % thermal noise psd
rx.Id = 10e-9; % dark current
rx.R = 1; % responsivity
% Electric Lowpass Filter
% rx.elefilt = design_filter('bessel', 5, 0.5*mpam.Rs/(sim.fs/2));
rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);
% Optical Bandpass Filter
rx.optfilt = design_filter('fbg', 0, sim.M*rx.elefilt.fcnorm);

% KLSE Fourier Series Expansion (done here because depends only on filters
% frequency response)
[rx.U_fourier, rx.D_fourier, rx.Fmax_fourier] = klse_fourier(rx, sim, sim.Mct*(mpam.M^sim.L + 2*sim.L)); 

GainsdB = 5:2.5:20;

figure, hold on, grid on
legends = {};
for k= 1:length(GainsdB)

    %% SOA
    % soa(GaindB, NF, lambda, maxGaindB)
    soaG = soa(GainsdB(k), 9, 1310e-9, 20); 

    % BER
    ber_soa = soa_ber(mpam, tx, soaG, rx, sim);
    
%     plot(tx.PtxdBm, log10(ber_soa.count), '-o')
    plot(tx.PtxdBm, log10(ber_soa.est), '-')
%     plot(tx.PtxdBm, log10(ber_soa.gauss), '--')
    legends = [legends, sprintf('Gain = %.1f dB', GainsdB(k))];
end
xlabel('Received Power (dBm)')
ylabel('log(BER)')
legend(legends{:})





    