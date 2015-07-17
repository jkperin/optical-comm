%% Vary Gain of Optical Amplifier
clear, clc

addpath ../f % general functions
addpath f

% Simulation parameters
sim.Nsymb = 2^14; % Number of symbols in montecarlo simulation
sim.Mct = 17;     % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
sim.L = 2;        % de Bruijin sub-sequence length (ISI symbol length)
sim.Me = 16; % Number of used eigenvalues
sim.BERtarget = 1e-4; 
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

sim.shot = false; % include shot noise in montecarlo simulation 
sim.RIN = false; % include RIN noise in montecarlo simulation
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
if mpam.M == 8
    tx.PtxdBm = -22:2:-8;
elseif mpam.M == 4
    tx.PtxdBm = -26:1:-10;
elseif mpam.M == 16
    tx.PtxdBm = -18:2:-2;
end

tx.lamb = 1310e-9; % wavelength
tx.alpha = 2; % chirp parameter
tx.RIN = -150;  % dB/Hz
tx.rexdB = -Inf;  % extinction ratio in dB. Defined as Pmin/Pmax

% Modulator frequency response
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
soa_opt = soa(20, 9, tx.lamb, Inf);
% soa_opt.optimize_gain(mpam, tx, fiber, rx, sim)

GainsdB = [5:19 soa_opt.GaindB];

soaG = soa(10, 9, 1310e-9, 20); 

figure, hold on, grid on
legends = {};
for k= 1:length(GainsdB)
    % Selected gain
    soaG.GaindB = GainsdB(k);

    %% Equally-spaced levels
    mpam.level_spacing = 'equally-spaced';
        
    ber.eq_spaced = soa_ber(mpam, tx, fiber, soaG, rx, sim);
    
    Preq.eq_spaced(k) = interp1(log10(ber.eq_spaced.est), tx.PtxdBm, log10(sim.BERtarget), 'spline');
    
    
    %% Optimized levels
    mpam.level_spacing = 'optimized';
        
    ber.optimized = soa_ber(mpam, tx, fiber, soaG, rx, sim);
    
    Preq.optimized(k) = interp1(log10(ber.optimized.est), tx.PtxdBm, log10(sim.BERtarget), 'spline');
    
%     plot(tx.PtxdBm, log10(ber_soa.count), '-o')
    plot(tx.PtxdBm, log10(ber.eq_spaced.est), '-')
    plot(tx.PtxdBm, log10(ber.optimized.est), '--')
    legends = [legends, sprintf('Gain = %.1f dB (equally-spaced)', GainsdB(k))];
    legends = [legends, sprintf('Gain = %.1f dB (optimized)', GainsdB(k))];
end
xlabel('Received Power (dBm)')
ylabel('log(BER)')
legend(legends{:})
axis([tx.PtxdBm([1 end]) -8 0])

figure(2), box on, grid on, hold on
plot(GainsdB, Preq.eq_spaced, '-o');
xlabel('SOA Gain (dB)')
ylabel('Required Transmitted Power (dBm)')

figure(3), box on, grid on, hold on
plot(GainsdB, Preq.optimized, '-o');
xlabel('SOA Gain (dB)')
ylabel('Required Transmitted Power (dBm)')



    