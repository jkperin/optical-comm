function margin_vs_gain()
clc, close all

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
tx.PtxdBm = -25:1:-10;

tx.lamb = 1310e-9; % wavelength
tx.alpha = 0; % chirp parameter
tx.RIN = -140;  % dB/Hz
tx.rexdB = -10;  % extinction ratio in dB. Defined as Pmin/Pmax

% Modulator frequency response
% tx.modulator.fc = 2*mpam.Rs; % modulator cut off frequency
% tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
% tx.modulator.h = @(t) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0));
% tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

%% b2b
b2b = fiber();

%% Receiver
rx.N0 = (30e-12).^2; % thermal noise psd
% Electric Lowpass Filter
% rx.elefilt = design_filter('bessel', 5, 1.25*mpam.Rs/(sim.fs/2));
rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);


Gains = 1:20;
GainsdB = 10*log10(Gains);

mpam_eqs = mpam;
mpam_opt = mpam;
mpam_opt.level_spacing = 'optimized';

ka = [0.1 0.25 0.5];

MargindB_eqs = zeros(length(ka), length(GainsdB));
Gopt_margin_eqs = zeros(length(ka), 1);
MargindB_opt = zeros(length(ka), length(GainsdB));
Gopt_margin_opt = zeros(length(ka), 1);

figure, hold on, grid on
for n = 1:length(ka)
    [MargindB_eqs(n, :), Gopt_margin_eqs(n)] = iterate(GainsdB, mpam_eqs, ka(n));
    [MargindB_opt(n, :), Gopt_margin_opt(n)] = iterate(GainsdB, mpam_opt, ka(n));
end
    
leg = {};
for n = 1:length(ka)
    hline(n) = plot(Gains, MargindB_eqs(n, :));
    leg = [leg sprintf('ka = %.2f', ka(n))];
end    

for n = 1:length(ka)
    plot(Gopt_margin_eqs(n), interp1(Gains, MargindB_eqs(n, :), Gopt_margin_eqs(n)), 'o', 'Color', get(hline(n), 'Color'));
    plot(Gains, MargindB_opt(n, :), '--', 'Color', get(hline(n), 'Color'))
    plot(Gopt_margin_opt(n), interp1(Gains, MargindB_opt(n, :), Gopt_margin_opt(n)), 'o', 'Color', get(hline(n), 'Color'));
end  

xlabel('APD Gain (Linear Units)')
ylabel(sprintf('Transmitted Optical Power (dBm) @ BER = %g', sim.BERtarget))
legend(leg)

% axis([Gains(1) Gains(end) -21 -19]);

function [MargindB, Gopt_margin] = iterate(GainsdB, mpam, ka)
    
    %% APD 
    % (GaindB, ka, GainBW, R, Id) 
    apdG = apd(10, 0.1, Inf, 1, 10e-9);
    pin = apd(0, 0, Inf, 1, 10e-9);
    
    PtxdBm_BERtarget = zeros(size(GainsdB));
    figure, hold on, grid on
    legends = {};
    for k= 1:length(GainsdB)

        apdG.GaindB = GainsdB(k);

        % BER
        ber_apd = apd_ber(mpam, tx, b2b, apdG, rx, sim);

        % Calculate power at the target BER
        PtxdBm_BERtarget(k) = interp1(log10(ber_apd.gauss), tx.PtxdBm, log10(sim.BERtarget));
    end
    
    ber_pin = apd_ber(mpam, tx, b2b, pin, rx, sim);
    
    
    PtxdBm_pin_BERtarget = interp1(log10(ber_pin.gauss), tx.PtxdBm, log10(sim.BERtarget));
    
    MargindB = PtxdBm_BERtarget - PtxdBm_pin_BERtarget;

      
    % Find optimal margin
    Gopt_margin = apdG.optGain(mpam, tx, b2b, rx, sim, 'margin'); 

    end
end

    
    