%% BER as a function of the clipping ratio at the transmitter only. 
%% Quantization is not included. If quantization is included, then clipping
%% ratio at the receiver must also be provided.
clear, clc, close all
format compact

addpath f           % functions path

%% Parameters to swipe
Fnl = 20e9;

RCLIP = 3:0.25:5;

%% Simulation parameters
sim.type = 'palloc';                % Type of compensation at the transmitter {'preemphasis', 'palloc'}
sim.BERtarget = 1.8e-4;                    % Target BER
sim.Nsymb = 2^14;                   % number of OFDM symbols
sim.Mct = 8;                        % oversampling rate for emulation of continuous time 
                                    % Note: Mct has to be at least 4 so DT filter design approximates CT well enough
sim.Navg = 2;                       % Number of noise realizations 
sim.verbose = false;                % show intermediate results
                                    
sim.shot = false;
sim.RIN = false;
sim.quantiz = false;                % include quantization at both transmitter and receiver
sim.ENOB = 5;                       % effective number of bits

%% OFDM
% Number of subcarriers, Number of used subcarriers, Constellation size,
% bit rate
ofdm = ofdm(64, 52, 16, 112e9);

%% Transmitter parameters 
tx.lamb = 1296.2e-9;
tx.kappa = 1;    % current to optical power conversion (dc slope)
tx.alpha = 2;    % modulator chirp
tx.RIN = -150;   % RIN in dB/Hz. Only used if sim.RIN is true
tx.rexdB = -Inf;   % exctinction ratio in dB. Defined as Pmin/Pmax

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
fiber = fiber(); % back-to-back
                                             
%% Receiver parameters
rx.R = 1;                           % responsivity
rx.NEP = 30e-12;                    % Noise equivalent power of the TIA at the receiver (A/sqrt(Hz))
rx.Sth = 2*rx.R^2*rx.NEP^2/2;       % two-sided psd of thermal noise at the receiver (Sth = N0/2)

% Antialiasing filter
rx.filter = design_filter('gaussian', 4, 1/(ofdm.Ms*sim.Mct));             
    
%% Initiliaze variables
PtxdBm = zeros(size(RCLIP));
bercount = zeros(size(RCLIP));
berest = zeros(size(RCLIP));

rng_seed = rng;
for k = 1:length(RCLIP)
    tx.rclip = RCLIP(k);                    % Clipping ratio of ACO-OFDM
     
    PtxdBmv = zeros(1, sim.Navg);
    berc = zeros(1, sim.Navg);
    bere = zeros(1, sim.Navg);                               
    
    %%
    rng(rng_seed);
    for kk = 1:sim.Navg
        [ber, Ptx] = dc_ofdm(ofdm, tx, fiber, rx, sim);
        
        PtxdBmv(kk) = 10*log10(Ptx/1e-3);
        berc(kk) = log10(ber.count);
        bere(kk) = log10(ber.est);
    end
   
    % format some important variables
    PtxdBm(k) = mean(PtxdBmv);
    bercount(k) = mean(berc);
    berest(k) = mean(bere); 
end

%% Required power for DC-OFDM AWGN (with ideal filters, without CP, without dc bias)
snr = fzero(@(x) berqam(ofdm.CS, x) - sim.BERtarget, 20);      % snr to achieve target BER;
snrl = 10^(snr/10);

Pnrx = snrl*ofdm.Ms*ofdm.Rs*rx.Sth/ofdm.Nc;
Pnawgn = (1-qfunc(tx.rclip))^2*Pnrx./abs(rx.R*tx.kappa)^2*ones(1, ofdm.Nu/2);
Pawgn = tx.kappa*sqrt(2*sum(Pnawgn))*(tx.rclip*(1 - qfunc(tx.rclip)) + 1/sqrt(2*pi)*exp(-tx.rclip^2/2));
PawgndBm = 10*log10(Pawgn/1e-3);

%% Required power for OOK AWGN
PookdBm = 10*log10(1/rx.R*qfuncinv(sim.BERtarget)*sqrt(rx.Sth*ofdm.Rb)/1e-3);

% Same as solving for square constellations
% qfuncinv((sim.BERtarget*log2(ofdm.CS)/4)/(1-1/sqrt(ofdm.CS)))/sqrt(3*pi*log2(ofdm.CS)/((ofdm.CS-1)*4*rx.Sth*ofdm.Rb))

% Power penalty with respect to NRZ-OOK assuming no bandwidth limitation.
% Oversampling penalty is included
PP_ofdm_dc = @ (M, r) 10*log10(r*sqrt(2*(M-1)./(3*log2(M)))); 

power_pen_ook = PtxdBm - PookdBm;

%% Figures
figure
plot(RCLIP, bercount, '-sk', RCLIP, berest, '-or')
hold on
plot(RCLIP, log10(sim.BERtarget)*ones(size(RCLIP)), 'k')
% plot(RCLIP, berint.', 'xk', 'MarkerSize', 6)
xlabel('Clipping ratio', 'FontSize', 12)
ylabel('log_{10}(BER)', 'FontSize', 12)
legend('BER counted', 'BER estimaded', 'Target BER')

figure
plot(RCLIP, power_pen_ook, '-sk', 'LineWidth', 2)
hold on
plot(RCLIP, PP_ofdm_dc(ofdm.CS, tx.rclip)*ones(size(Fnl)), '-r', 'LineWidth', 2)
xlabel('Clipping ratio', 'FontSize', 12)
ylabel('Optical Power Penalty (dB)', 'FontSize', 12)
legend('Power penalty OFDM', 'Power penalty OFDM in AWGN')

