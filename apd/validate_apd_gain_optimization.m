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
tx.PtxdBm = -25:0.5:-10;

tx.lamb = 1310e-9; % wavelength
tx.alpha = 0; % chirp parameter
tx.RIN = -140;  % dB/Hz
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
% (GaindB, ka, GainBW, R, Id)  
apd_opt = apd(10.0851, 0.09, Inf, 1, 10e-9); % infinite gain x BW product
apd_opt.optimize_gain(mpam, tx, fiber, rx, sim);
        
%
GainsdB = sort([(6:0.5:13) apd_opt.GaindB]);
% GainsdB = apd_opt.GaindB;
Gains = 10.^(GainsdB/10);

% APD
apdG = apd(10, apd_opt.ka, Inf, 1, 10e-9);

% 
PtxdBm_BERtarget = zeros(size(GainsdB));
figure, hold on, grid on
legends = {};
for k= 1:length(GainsdB)

    apdG.GaindB = GainsdB(k);

    % BER
    ber_apd = apd_ber(mpam, tx, fiber, apdG, rx, sim);
     
    % Calculate power at the target BER
    PtxdBm_BERtarget(k) = interp1(log10(ber_apd.gauss), tx.PtxdBm, log10(sim.BERtarget));
    
%     plot(tx.PtxdBm, log10(ber_apd.count), '-o')
    plot(tx.PtxdBm, log10(ber_apd.gauss), '-')

    legends = [legends, sprintf('Gain = %.1f dB', GainsdB(k))];
end


%% Noise calculations
% Thermal noise
Deltaf = rx.elefilt.noisebw(sim.fs)/2; % electric filter one-sided noise bandwidth
varTherm = rx.N0*Deltaf; % variance of thermal noise

% RIN
if isfield(sim, 'RIN') && sim.RIN
    varRIN =  @(Plevel) 10^(tx.RIN/10)*Plevel.^2*Deltaf;
else
    varRIN = @(Plevel) 0;
end

Ptx = 1e-3*10.^(tx.PtxdBm/10);
rex = 10^(-abs(tx.rexdB)/10);
for k = 1:length(Ptx)
    mpam.adjust_levels(Ptx(k), tx.rexdB);
    
    Gapd(k) = apdG.optGain_analytical(mpam, rx.N0);
    
    apdG.Gain = Gapd(k);
    
    mpam.adjust_levels(apdG.Gain*Ptx(k), tx.rexdB);
    
    noise_std = @(Plevel) sqrt(varTherm + varRIN(Plevel) + apdG.varShot(Plevel/apdG.Gain, Deltaf));
    
    ber2(k) = mpam.ber_awgn(noise_std);
end
    

plot(tx.PtxdBm, log10(ber2), '-k')
    
    
    
    

xlabel('Received Power (dBm)')
ylabel('log(BER)')
legend(legends{:})
axis([tx.PtxdBm(1) tx.PtxdBm(end) -8 0])

figure, hold on, grid on
plot(Gains, PtxdBm_BERtarget)
plot(apd_opt.Gain, PtxdBm_BERtarget(GainsdB == apd_opt.GaindB), 'ok')
xlabel('APD Gain (Linear Units)')
ylabel(sprintf('Transmitted Optical Power (dBm) @ BER = %g', sim.BERtarget))
% axis([Gains(1) Gains(end) -21 -19]);
    