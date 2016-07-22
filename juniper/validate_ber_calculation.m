%% Validate BER calculation of pre-amplified system

clear, clc, close all

addpath f/ % Juniper project specific functions
addpath ../mpam % PAM
addpath ../f % general functions
addpath ../soa % for pre-amplifier 
addpath ../apd % for PIN photodetectors

%% Transmit power swipe
Tx.PtxdBm = -30:-20; % transmitter power range
% Tx.PtxdBm = -10; % transmitter power range
Rx.PrxdBm = 5; % constant power at the receiver

%% Simulation parameters
sim.Rb = 56e9;    % bit rate in bits/sec
sim.Nsymb = 2^10; % Number of symbols in montecarlo simulation
% For DACless simulation must make Tx.dsp.ros = sim.Mct and DAC.resolution = Inf
sim.Mct = 1;      % Oversampling ratio to simulate continuous time. Must be integer multiple of sim.ros.txDSP and numerator of sim.ros.rxDSP
sim.BERtarget = 1e-4; 
sim.Ndiscard = 128; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation
sim.Modulator = 'MZM'; % 'MZM' only

%% Pulse shape
pulse_shape = select_pulse_shape('rect', sim.Mct);
% pulse_shape = select_pulse_shape('rrc', sim.ros.txDSP, 0.2, 4);
Tx.pulse_shape = pulse_shape;

%% M-PAM
% PAM(M, bit rate, leve spacing : {'equally-spaced', 'optimized'}, pulse
% shape: struct containing properties of pulse shape 
mpam = PAM(4, sim.Rb, 'equally-spaced', pulse_shape);

%% Amplifier
% Class SOA characterizes amplifier in terms of gain and noise figure
% soa(GaindB: amplifier gain in dB, NFdB: noise figure in dB, lambda: operationa wavelength, maxGaindB: maximum amplifier gain) 
EDFA = soa(20, 5, 1550e-9, 20); 

dataTX = randi([0 mpam.M-1], [1 sim.Nsymb]);

%% AWGN approximation
% Note: this calculation includes noise enhancement due to equalization,
% but dominant noise in pre-amplified system is the signal-spontaneous beat
% noise, which is not Gaussian distributed
noiseBW = mpam.Rs/2;

% Noise std for intensity level Plevel
Npol = 2; % number of polarizations. Npol = 1, if polarizer is present, Npol = 2 otherwise.


%% 
mpamRef = mpam;
for k = 1:length(Tx.PtxdBm)
    mpam = mpam.adjust_levels(dBm2Watt(Tx.PtxdBm(k)), -Inf);
    
    x = mpam.signal(dataTX);
    
    [Eamp, OSNRdB(k)] = EDFA.amp(sqrt(x), mpam.Rs);
    
    Att = dBm2Watt(Rx.PrxdBm)/dBm2Watt(power_meter(Eamp));
    Eamp = Eamp*sqrt(Att);
    power_meter(Eamp)
    
    y = abs(Eamp).^2;
    
    y = y/(EDFA.Gain*Att);
    
    dataRX = mpam.demod(y(1, :));
    
    [~, ber_count(k)] = biterr(dataTX, dataRX);
   
    % AWGN approximation
    mpamRef = mpamRef.adjust_levels(dBm2Watt(Rx.PrxdBm), -Inf);
    noise_std = @(Plevel) sqrt(2*Att*Plevel*EDFA.N0*noiseBW);

    ber_awgn(k) = mpamRef.berAWGN(noise_std);
end

figure(1), box on, hold on
plot(OSNRdB, log10(ber_count));
plot(OSNRdB, log10(ber_awgn));
OSNRdBexp = [23,27,29,31,33,35,38,40,42,44,46,48];
BERcountexp = [0.00214091492392766,0.000493871972323897,0.000158199639915134,5.21978507334200e-05,1.52578332913074e-05,4.81826314462338e-06,2.40913157231169e-06,8.01344977410085e-07,8.01344977410085e-07,0,0,0];
plot(OSNRdBexp, log10(BERcountexp), '-ok', 'linewidth', 2)
legend('Counted', 'Gaussian approximation', 'Experiment')
axis([min(OSNRdB(1), OSNRdBexp(1)) OSNRdBexp(end) -8 0])


    
    



