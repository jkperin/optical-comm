%% BER vs transmitted power
clear, close all

format compact

addpath f           % functions path

% rng('shuffle');     % Reinitialize the random number generator used by rand, randi, and randn with a seed based on the current time

%% Parameters to change
sim.ENOB = 6;                       % Effective number of bits for DAC and ADC (only used if sim.quantiz is true)
sim.quantiz = ~true;                 % include quantization at both transmitter and receiver
sim.shot = ~true;                    % Include shot noise?
sim.RIN = true;                     % Include intensity noise?
sim.verbose = false;                    % Show all plots? It'll slow donw simulation
sim.Navg = 1;                     % Number of noise realizations
                      
sim.type = 'preemphasis';           % type of power allocation
sim.Pb = 1.8e-4;                    % Target BER
sim.Nsymb = 2^10;                   % number of OFDM symbols
sim.Mct = 4;                       % oversampling rate for emulation of continuous time 
                                    % Note: Mct has to be at least 4 so DT filter design approximates CT well enough./
                                                                      
%% Modulation parameters 
%% Modulation parameters 
% Number of subcarriers, Number of used subcarriers, Constellation size,
% bit rate
ofdm = ofdm(64, 52, 16, 112e9);

%% Fiber
fiber = fiber();

%% Transmitter parameters 
tx.lamb = 1296.2e-9;
tx.modulator.fc = 20e9; % modulator cut off frequency
tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
tx.modulator.h = @(t) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0));
tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

tx.kappa = 1;    % current to optical power conversion (dc slope)

tx.alpha = 1;    % modulator chirp

tx.rexdB = -Inf;   % exctinction ratio in dB. Defined as Pmin/Pmax
tx.RIN = -160;   % RIN in dB/Hz. Only used if sim.RIN is true

% Transmitter filter (ZOH + some smoothing filter)
tx.filter = design_filter('bessel', 5, 1/(ofdm.Ms*sim.Mct));

% Convolve with ZOH
bzoh = ones(1, sim.Mct)/sim.Mct;
tx.filter.num = conv(tx.filter.num, bzoh);
tx.filter.grpdelay = grpdelay(tx.filter.num, tx.filter.den, 1);
tx.filter.H = @(f) freqz(tx.filter.num, tx.filter.den, 2*pi*f).*exp(1j*2*pi*f*tx.filter.grpdelay);
tx.filter.noisebw = @(fs) noisebw(tx.filter.num, tx.filter.den, 2^15, fs); % equivalent two-sided noise bandwidth over larger number of samples (2^15) given a sampling frequency fs 

%% Receiver parameters
rx.R = 1;                           % responsivity
rx.NEP = 20e-12;                    % Noise equivalent power of the TIA at the receiver (A/sqrt(Hz))
rx.Sth = rx.R^2*rx.NEP^2/2;       % two-sided psd of thermal noise at the receiver (Sth = N0/2)

% Antialiasing filter
rx.filter = design_filter('gaussian', 4, 1/(ofdm.Ms*sim.Mct));        

%% define clipping ratios based on optimization results
% this script creates the variables sim.rcliptx, sim.rcliprx, RCLIPTX and
% RCLIPRX based on optimized values and on whether or not quantization is
% on
select_clipping_ratio
tx.rclip = sim.rcliptx;
rx.rclip = sim.rcliptx;

PrxdBm = -10:-4;
Prx = 1e-3*10.^(PrxdBm/10);

% Initiliaze variables
bercount1 = zeros(size(Prx));
berest1 = zeros(size(Prx));

bercount2 = zeros(size(Prx));
berest2 = zeros(size(Prx));

%% Simulation 1: preemphasis, B2B
sim1 = sim;
sim1.type = 'preemphasis';

%% Simulation 2: palloc, B2B
sim2 = sim;
sim2.type = 'palloc';

clear sim
for k = 1:length(Prx)
    
    tx.Ptx = Prx(k);
    
    for kk = 1:sim1.Navg
        % generate, propagate, and detect
        ber1 = dc_ofdm(ofdm, tx, fiber, rx, sim1);
        
        ber2 = dc_ofdm(ofdm, tx, fiber, rx, sim2);

        berc1(kk) = ber1.count;
%         bere1(kk) = ber1.est;
        
        berc2(kk) = ber2.count;
%         bere2(kk) = ber2.est;
    end
    
    bercount1(k) = log10(mean(berc1));
    berest1(k) = log10(mean(berc1));
    
    bercount2(k) = log10(mean(berc2));
    berest2(k) = log10(mean(berc2));
        
end

%% "Experimental study of PAM-4, CAP-16, and DMT for 100 Gb/s Short Reach Optical Transmission Systems"
oepaper.PrxdBm = -10:-5;
oepaper.BER = [9e-2, 5e-2, 2e-2, 4.5e-3, 1e-3, 7e-5];
oepaper.PrxdBmexp = -10:-1;
oepaper.BERexp = [4e-2, 2e-2, 1.3e-2, 8e-3, 6e-3, 4e-3, 3e-3, 3e-3, 3e-3, 3e-3];

%% Figures
figure, hold on, grid on, box on
plot(PrxdBm, berest1, '-ob')
plot(PrxdBm, bercount1, '--xb')
plot(PrxdBm, berest2, '-or')
plot(PrxdBm, bercount2, '--xr')
plot(oepaper.PrxdBm, log10(oepaper.BER), '-sm')
plot(oepaper.PrxdBmexp, log10(oepaper.BERexp), '-ok')
xlabel('P_{rx} (dBm)')
ylabel('log_{10}(BER)')
legend('Preemphasis (Analysis)', 'Preemphasis (Montecarlo)', ...
    'Palloc (Analysis)', 'Palloc (Montecarlo)', 'oe-23-2-1176 Simulation B2B',...
    'oe-23-2-1176 Experiment B2B');
axis([-10 -2 -5 0])
% saveas(gca, sprintf('comparison_with_oepaper_%dkm', fiber.L/1e3), 'png')