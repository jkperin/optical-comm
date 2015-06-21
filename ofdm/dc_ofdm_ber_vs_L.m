%% BER vs fiber length
clear, clc, close all
format compact

addpath f           % functions path

%% Parameters to change
L = (0:1:5)*1e3;                   % fiber length

Fnl = 30e9;                         % modulator bandwidth

tx.channel = -3;                    % channel number. 0 corresponds to zero-wavelength dispersion.

tx.dlamb = 20e-9;                   % WDM channel spacing (m)   

%% Simulation parameters
sim.type = 'palloc';                % Type of compensation at the transmitter {'preemphasis', 'palloc'}
sim.BERtarget = 1.8e-4;                    % Target BER
sim.Nsymb = 2^12;                   % number of OFDM symbols
sim.Mct = 8;                        % oversampling rate for emulation of continuous time 
                                    % Note: Mct has to be at least 4 so DT filter design approximates CT well enough
sim.Navg = 1;                       % Number of noise realizations 
sim.verbose = false;                % show intermediate results
                                    
sim.shot = true;
sim.RIN = true;
sim.quantiz = true;                % include quantization at both transmitter and receiver
sim.ENOB = 6;                       % effective number of bits

%% OFDM
% Number of subcarriers, Number of used subcarriers, Constellation size,
% bit rate
ofdm = ofdm(64, 52, 16, 112e9);

%% Transmitter parameters 
tx.lamb = 1310e-9 + tx.channel*tx.dlamb;
tx.kappa = 1;    % current to optical power conversion (dc slope)
tx.alpha = 2;    % modulator chirp
tx.RIN = -150;   % RIN in dB/Hz. Only used if sim.RIN is true
tx.rexdB = -10;   % exctinction ratio in dB. Defined as Pmin/Pmax

% Modulator frequency response
tx.modulator.fc = Fnl; % modulator cut off frequency
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
att = @(x) 0.35; % constant attenuation at all wavelengths
fiber = fiber(0, att); % smf28      

%% Receiver parameters
rx.R = 1;                           % responsivity
rx.NEP = 30e-12;                    % Noise equivalent power of the TIA at the receiver (A/sqrt(Hz))
rx.Sth = 2*rx.R^2*rx.NEP^2/2;       % two-sided psd of thermal noise at the receiver (Sth = N0/2)

% Antialiasing filter
rx.filter = design_filter('gaussian', 4, 1/(ofdm.Ms*sim.Mct));       

%% define clipping ratios based on optimization results
[tx.rclip, rx.rclip] = select_clipping_ratio(sim.type, ofdm.CS, tx.modulator.fc);

% Initiliaze variables
PtxdBm = zeros(size(L));
PtxedBm = zeros(size(L));
bercount = zeros(size(L));
berest = zeros(size(L));
PookdBm = zeros(size(L));

for k = 1:length(L)
    
    fiber.L = L(k);                    
    
    for kk = 1:sim.Navg
        % generate, propagate, and detect
        [ber, Ptx(kk), Ptx_est(kk)] = dc_ofdm(ofdm, tx, fiber, rx, sim);

        % format some important variables
        berc(kk) = ber.count;
        bere(kk) = ber.est;
    end
    
    PtxdBm(k) = 10*log10(mean(Ptx)/1e-3);
    PtxedBm(k) = 10*log10(mean(Ptx_est)/1e-3);
    bercount(k) = log10(mean(berc));
    berest(k) = log10(mean(bere));
    
    %% Required power for OOK AWGN
    PookdBm(k) = 10*log10(1/rx.R*qfuncinv(sim.BERtarget)*sqrt(rx.Sth*ofdm.Rb)/1e-3) + fiber.att(tx.lamb)*fiber.L/1e3;

    % Same as solving for square constellations
    % qfuncinv((sim.BERtarget*log2(ofdm.CS)/4)/(1-1/sqrt(ofdm.CS)))/sqrt(3*pi*log2(ofdm.CS)/((ofdm.CS-1)*4*rx.Sth*ofdm.Rb))
    
end

power_pen_ook_m = PtxdBm - PookdBm;
power_pen_ook_e = PtxedBm - PookdBm;


%% Figures
figure, hold on, grid on
plot(L/1e3, bercount, '-sk', L/1e3, berest, '-or')
plot(L/1e3, log10(sim.BERtarget)*ones(size(L)), 'k')
xlabel('Length (km)', 'FontSize', 12)
ylabel('log_{10}(BER)', 'FontSize', 12)
legend('BER counted', 'BER estimaded', 'Target BER', 'Location', 'NorthWest')
axis([L([1 end])/1e3 -4.1 -3])

figure, hold on, grid on
plot(L/1e3, power_pen_ook_m, '-sk', 'LineWidth', 2)
plot(L/1e3, power_pen_ook_e, '-xb', 'LineWidth', 2)
xlabel('Length (km)', 'FontSize', 12)
ylabel('Optical Power Penalty (dB)', 'FontSize', 12)
legend('Measured', 'Estimated', 'Location', 'NorthEast')
axis([L([1 end])/1e3 0 18])
