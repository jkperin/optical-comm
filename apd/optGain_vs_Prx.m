%% Validate APD Gain Optimization
clear, clc, close all

addpath ../mpam
addpath ../f
addpath f

% Simulation parameters
sim.Nsymb = 2^14; % Number of symbols in montecarlo simulation
sim.Mct = 9;    % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
sim.L = 2;        % de Bruijin sub-sequence length (ISI symbol length)
sim.BERtarget = 1e-4; 
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation
sim.awgn = true;
%
sim.shot = true; % include shot noise. Only included in montecarlo simulation (except for APD)
sim.RIN = true; % include RIN noise. Only included in montecarlo simulation
sim.verbose = ~true; % show stuff

% M-PAM
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
tx.PtxdBm = -20:1:0;

tx.lamb = 1310e-9; % wavelength
tx.alpha = 0; % chirp parameter
tx.RIN = -145;  % dB/Hz
tx.rexdB = -10;  % extinction ratio in dB. Defined as Pmin/Pmax

% Modulator frequency response
% tx.modulator.fc = 2*mpam.Rs; % modulator cut off frequency
% tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
% tx.modulator.h = @(t) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0));
% tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

%% Fiber
fiber = fiber();

%% Receiver
rx.N0 = (30e-12).^2; % thermal noise psd
% Electric Lowpass Filter
% rx.elefilt = design_filter('bessel', 5, 1.25*mpam.Rs/(sim.fs/2));
rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);

%% APD 
apdG = apd(10, 0.1, Inf, 1, 10e-9);

% Optimize APD gain: Given a certain input power calculates 
% the APD gain that leads to the minimum BER
Gopt = apdG.optGain(mpam, tx, fiber, rx, sim);

% Optimizes APD gain assuming equally-spaced levels and that highest level
% dominates BER
Gopt_analytical = zeros(size(tx.PtxdBm));
for k = 1:length(tx.PtxdBm)
    Ptx = 1e-3*10^(tx.PtxdBm(k)/10);
    
    mpam = mpam.adjust_levels(Ptx, tx.rexdB);
    
    Gopt_analytical(k) = apdG.optGain_analytical(mpam, rx.N0);    
end

% Optimize APD gain
% 1) For a given APD gain, optimizes level spacing for the target BER
% 2) Adjust power to the given received power
% 3) Calculate BER
% 4) Update APD gain until minimum BER is reached
mpam.level_spacing = 'optimized'
Gopt_opt = apdG.optGain(mpam, tx, fiber, rx, sim);
Gopt_opt(:) = 10;

figure, hold on
for kk = -2:2:2
    for k = 1:length(Gopt_opt)
        apdG.Gain = Gopt_opt(k) + kk;

        tx.Ptx = 1e-3*10^(tx.PtxdBm(k)/10);

        noise_std = apdG.stdNoise(rx.elefilt.noisebw(sim.fs)/2, rx.N0, tx.RIN);

        mpam = mpam.optimize_level_spacing_gauss_approx(sim.BERtarget, tx.rexdB, noise_std, sim.verbose);

        % Adjust levels to power after the APD, which is required by noise_std
%         link_gain = apdG.R*apdG.Gain*fiber.link_attenuation(tx.lamb);

%         mpam = mpam.adjust_levels(tx.Ptx*link_gain, tx.rexdB);

%         bergauss(k) = mpam.ber_awgn(noise_std);
        
        [~, bergauss(k)] = ber_apd_doubly_stochastic(mpam, tx, fiber, apdG, rx, sim);
    end
    
    plot(tx.PtxdBm, log10(bergauss))
end

xlabel('Received Power (dBm)')
ylabel('log(BER)')
legend('Gopt-2', 'Gopt', 'Gopt+2')
axis([tx.PtxdBm(1) tx.PtxdBm(end) -8 0])

figure, hold on, box on, grid on
plot(tx.PtxdBm, Gopt)
plot(tx.PtxdBm, Gopt_analytical, '--')
plot(tx.PtxdBm, Gopt_opt)
legend('Gopt equally-spaced', 'Gopt analytical equally-spaced', 'Gopt optimized level spacing')
xlabel('Received Power (dBm)')
ylabel('Gain')
    
    