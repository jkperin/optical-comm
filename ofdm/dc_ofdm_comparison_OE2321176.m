%% DC-OFDM BER x transmitted power
% Results are compared with the results presented in OE-23-2-1176 paper
clear, clc, close all
format compact

addpath f           % functions path

%% Transmitted power to swipe
PrxdBm = -10:-4;
Prx = 1e-3*10.^(PrxdBm/10);

%% Simulation parameters
sim.ENOB = 6;                       % Effective number of bits for DAC and ADC (only used if sim.quantiz is true)
sim.quantiz = false;                 % include quantization at both transmitter and receiver
sim.shot = false;                    % Include shot noise?
sim.RIN = true;                     % Include intensity noise?
sim.verbose = false;                    % Show all plots? It'll slow donw simulation
sim.Navg = 1;                     % Number of noise realizations
                      
sim.BERtarget = 1.8e-4;            % Target BER
sim.Nsymb = 2^10;                   % number of OFDM symbols
sim.Mct = 4;                       % oversampling rate for emulation of continuous time 
                                    % Note: Mct has to be at least 4 so DT filter design approximates CT well enough./
                                                                      
%% OFDM
% Number of subcarriers, Number of used subcarriers, Constellation size,
% bit rate
ofdm = ofdm(64, 52, 16, 112e9);

%% Transmitter parameters 
tx.lamb = 1296.2e-9;
tx.kappa = 1;    % current to optical power conversion (dc slope)
tx.alpha = 2;    % modulator chirp
tx.RIN = -160;   % RIN in dB/Hz. Only used if sim.RIN is true
tx.rexdB = -Inf;   % exctinction ratio in dB. Defined as Pmin/Pmax

% Modulator frequency response
tx.modulator.fc = 20e9; % modulator cut off frequency
tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
tx.modulator.h = @(t) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0));
tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

% Transmitter filter (ZOH + some smoothing filter)
tx.filter = design_filter('bessel', 5, 1/(ofdm.Ms*sim.Mct));

% Convolve with ZOH
bzoh = ones(1, sim.Mct)/sim.Mct;
tx.filter.num = conv(tx.filter.num, bzoh);
tx.filter.grpdelay = grpdelay(tx.filter.num, tx.filter.den, 1);
tx.filter.H = @(f) freqz(tx.filter.num, tx.filter.den, 2*pi*f).*exp(1j*2*pi*f*tx.filter.grpdelay);
tx.filter.noisebw = @(fs) noisebw(tx.filter.num, tx.filter.den, 2^15, fs); % equivalent two-sided noise bandwidth over larger number of samples (2^15) given a sampling frequency fs 

%% Fiber
fiber = fiber(); % back-to-back

%% Receiver parameters
rx.R = 1;                           % responsivity
rx.NEP = 20e-12;                    % Noise equivalent power of the TIA at the receiver (A/sqrt(Hz))
rx.Sth = 2*rx.R^2*rx.NEP^2/2;       % two-sided psd of thermal noise at the receiver (Sth = N0/2)

% Antialiasing filter
rx.filter = design_filter('gaussian', 4, 1/(ofdm.Ms*sim.Mct));        

%% Initiliaze variables
bercount1 = zeros(size(Prx));
berest1 = zeros(size(Prx));

bercount2 = zeros(size(Prx));
berest2 = zeros(size(Prx));

%% Simulation 1: preemphasis, B2B
sim1 = sim;
sim1.type = 'preemphasis';
% [tx.rclip, rx.rclip] = select_clipping_ratio(sim1.type, ofdm.CS, tx.modulator.fc);

%% Simulation 2: palloc, B2B
sim2 = sim;
sim2.type = 'palloc';
[tx.rclip, rx.rclip] = select_clipping_ratio(sim2.type, ofdm.CS, tx.modulator.fc);

clear sim
for k = 1:length(Prx)
    
    tx.Ptx = Prx(k)*10^(fiber.att(tx.lamb)*fiber.L/1e4);
    
    for kk = 1:sim1.Navg
        % generate, propagate, and detect
        ber1 = dc_ofdm(ofdm, tx, fiber, rx, sim1);
        
        ber2 = dc_ofdm(ofdm, tx, fiber, rx, sim2);

        berc1(kk) = ber1.count;
        bere1(kk) = ber1.est;
        
        berc2(kk) = ber2.count;
        bere2(kk) = ber2.est;
    end
    
    bercount1(k) = log10(mean(berc1));
    berest1(k) = log10(mean(bere1));
    
    bercount2(k) = log10(mean(berc2));
    berest2(k) = log10(mean(bere2));  
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