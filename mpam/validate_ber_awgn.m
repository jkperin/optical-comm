%% Compare performance of M-PAM in 3 different scenarios
% 1. PIN receiver, no amplifier
% 2. SOA & PIN receiver
% 3. APD
clear, clc, close all

addpath ../f/

%% Simulation parameters
sim.Nsymb = 2^20; % Number of symbols in montecarlo simulation
sim.Mct = 15;     % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
sim.BERtarget = 1e-4; 
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

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

%% Receiver
rx.N0 = (20e-12).^2; % thermal noise psd

matchedfilt = design_filter('matched', mpam.pshape, 1/sim.Mct);

PtxdBm = -22:1:-10;
Ptx = 1e-3*10.^(PtxdBm/10);

ndiscard = [1:sim.Ndiscard sim.Nsymb-(0:sim.Ndiscard-1)];

for k = 1:length(Ptx)

    dataTX = randi([0 mpam.M-1], 1, sim.Nsymb);
    
    mpam.adjust_levels(Ptx(k), -Inf);

    x = mpam.mod(dataTX, sim.Mct);
    
    x([1:sim.Ndiscard*sim.Mct end-(0:sim.Ndiscard*sim.Mct-1 )]) = 0;

    p = x + sqrt(rx.N0*sim.fs/2)*randn(size(x));

    pf = real(ifft(fft(p).*ifftshift(matchedfilt.H(f/sim.fs))));

    y = pf(floot(sim.Mct/2):sim.Mct:end);
        
    dataRX = mpam.demod(y);
    
    dataTX(ndiscard) = [];
    dataRX(ndiscard) = [];
    [~, ber.count(k)] = biterr(dataTX, dataRX);
    
    ber.awgn(k) = mpam.ber_awgn(@(P) sqrt(rx.N0*matchedfilt.noisebw(sim.fs)));
    
end

figure, hold on, grid on, box on
plot(PtxdBm, log10(ber.count), '-o')
plot(PtxdBm, log10(ber.awgn), '-')
axis([PtxdBm([1 end]) -10 0])
