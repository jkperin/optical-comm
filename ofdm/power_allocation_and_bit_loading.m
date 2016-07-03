%% Power allocation for a given channel
clear, clc

addpath ../f/
addpath f/

%% Parameters to change
sim.ENOB = 6;                       % Effective number of bits for DAC and ADC (only used if sim.quantiz is true)
sim.quantiz = false;                 % include quantization at both transmitter and receiver
sim.shot = true;                    % Include shot noise?
sim.RIN = true;                     % Include intensity noise?
sim.verbose = false;                    % Show all plots? It'll slow donw simulation
sim.Navg = 1;                     % Number of noise realizations
                      
sim.type = 'preemphasis';           % type of power allocation
sim.BERtarget = 1.8e-4;                    % Target BER
sim.Nsymb = 2^10;                   % number of OFDM symbols
sim.Mct = 4;                       % oversampling rate for emulation of continuous time 
                                    % Note: Mct has to be at least 4 so DT filter design approximates CT well enough./
                                                                      
%% Modulation parameters 
%% Modulation parameters 
% Number of subcarriers, Number of used subcarriers, Constellation size,
% bit rate
ofdm = ofdm(512, 416, 16, 56e9); 

%% Transmitter parameters 
tx.modulator.fc = 30e9; % modulator cut off frequency
tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
tx.modulator.h = @(t) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0));
tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

tx.lamb = 1550e-9;

tx.kappa = 1;    % current to optical power conversion (dc slope)
tx.alpha = 0;    % modulator chirp

tx.rexdB = -Inf;   % exctinction ratio in dB. Defined as Pmin/Pmax
tx.RIN = -150;   % RIN in dB/Hz. Only used if sim.RIN is true
tx.rclip =6;
rx.rclip = 6;

% Transmitter filter (ZOH + some smoothing filter)
tx.filter = design_filter('bessel', 5, 1/(ofdm.Ms*sim.Mct));

% Convolve with ZOH
bzoh = ones(1, sim.Mct)/sim.Mct;
tx.filter.num = conv(tx.filter.num, bzoh);
tx.filter.grpdelay = grpdelay(tx.filter.num, tx.filter.den, 1);
tx.filter.H = @(f) freqz(tx.filter.num, tx.filter.den, 2*pi*f).*exp(1j*2*pi*f*tx.filter.grpdelay);
tx.filter.noisebw = @(fs) noisebw(tx.filter.num, tx.filter.den, 2^15, fs); % equivalent two-sided noise bandwidth over larger number of samples (2^15) given a sampling frequency fs


%% Amplifier
% Class SOA characterizes amplifier in terms of gain and noise figure
% soa(GaindB: amplifier gain in dB, NFdB: noise figure in dB, lambda: operationa wavelength, maxGaindB: maximum amplifier gain)
EDFA = soa(20, 5, 1550e-9, 20); 

%% Receiver parameters
rx.R = 1;                           % responsivity
rx.NEP = 20e-12;                    % Noise equivalent power of the TIA at the receiver (A/sqrt(Hz))
rx.Sth = rx.R^2*rx.NEP^2/2;       % two-sided psd of thermal noise at the receiver (Sth = N0/2)

% Antialiasing filter
rx.filter = design_filter('gaussian', 4, 1/(ofdm.Ms*sim.Mct));    

%% Photodiode
% pin(R: Responsivity (A/W), Id: Dark current (A), BW: bandwidth (Hz))
% PIN frequency response is modeled as a first-order system
Rx.PD = pin(1, 10e-9, Inf);

%% Fiber
Fiber = fiber(1e3);

%% Calculate cyclic prefix ##!! fiber has to be included here
ofdm.cyclic_prefix(@(f, fs) Fiber.Himdd(f, tx.lamb), true);

%% Time and frequency scales
sim.fs = sim.Mct*ofdm.fs;                                     % sampling frequency to emulate continuous time (Hz)  

sim.N = sim.Mct*(ofdm.Nc + ofdm.Npre_os)*sim.Nsymb;           % total number of points simulated in continuous time
         
dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

sim.t = t;
sim.f = f; 

PrxdBm = -10;
Gch = Fiber.Himdd(ofdm.fc, tx.lamb);
% Noise bandwidth
noiseBW = ofdm.fs/2;
BWopt = 50e9; % optical filter noise bandwidth; Not divided by 2 because optical filter is a bandpass filter

% Thermal noise
varTherm = rx.Sth*noiseBW; % variance of thermal noise

% Noise std for intensity level Plevel
Npol = 2; % number of polarizations. Npol = 1, if polarizer is present, Npol = 2 otherwise.
varNoise =  sqrt(varTherm + Rx.PD.varShot(dBm2Watt(PrxdBm), noiseBW)...
    + Rx.PD.R^2*EDFA.varAWGN(dBm2Watt(PrxdBm)/EDFA.Gain, noiseBW, BWopt, Npol));
% Note: Plevel is divided by amplifier gain to obtain power at the amplifier input

%% 
[Pn, CS] = palloc(ofdm, Gch, varNoise, sim.BERtarget);

figName =sprintf('Power allocation and bit loading for L = %d km of SMF', Fiber.L/1e3);
figure(1), clf
subplot(211), hold on, box on
stem(ofdm.fc/1e9, Pn/Pn(1))
plot(ofdm.fc/1e9, abs(Fiber.Himdd(ofdm.fc, tx.lamb, tx.alpha, 'large signal')).^2, 'LineWidth', 2)
xlabel('Subcarrier frequency (GHz)', 'FontSize', 14)
ylabel('Normalized allocated power', 'FontSize', 14)
title(figName)
set(gca, 'FontSize', 14)
grid on

subplot(212), hold on, box on
stem(ofdm.fc/1e9, CS)
xlabel('Subcarrier frequency (GHz)', 'FontSize', 14)
ylabel('Constellation size', 'FontSize', 14)
set(gca, 'FontSize', 14)
set(gca, 'ytick', [4 16 32 64])
grid on
saveas(gca, ['../juniper/figs/' figName], 'png')