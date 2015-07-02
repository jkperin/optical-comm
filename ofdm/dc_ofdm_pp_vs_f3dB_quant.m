%% Power penalty vs modulator cutoff frequency assuming quantization
clear, clc, close all
format compact

addpath f           % functions path

%% Parameters to change
Fnl = (20:5:50)*1e9;

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
tx.lamb = 1310e-9;
tx.kappa = 1;    % current to optical power conversion (dc slope)
tx.alpha = 2;    % modulator chirp
tx.RIN = -150;   % RIN in dB/Hz. Only used if sim.RIN is true
tx.rexdB = -10;   % exctinction ratio in dB. Defined as Pmin/Pmax

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
rx.Sth = rx.R^2*rx.NEP^2/2;       % two-sided psd of thermal noise at the receiver (Sth = N0/2)

% Antialiasing filter
rx.filter = design_filter('gaussian', 4, 1/(ofdm.Ms*sim.Mct));       

%% Iterate simulation varying the cut-off frequency of the 2nd-order filter of the transmitter
% Note: cut-off frequency changes the cyclic prefix length and thus chnages
% the sampling rate. 

% Initiliaze variables
PtxdBm = zeros(size(Fnl));
PtxedBm = zeros(size(Fnl));
bercount = zeros(size(Fnl));
berest = zeros(size(Fnl));
berint = zeros(length(Fnl), 2);

PP = zeros(length(Fnl), 2);

fprintf('------ %s, CS = %d, ENOB = %d ------\n', sim.type, ofdm.CS, sim.ENOB)

for k = 1:length(Fnl)
    fprintf('-------------- fnl = %d GHz --------------\n', Fnl(k)/1e9)   
    % Modulator frequency response
    tx.modulator.fc = Fnl(k); % modulator cut off frequency
    tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
    tx.modulator.h = @(t) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0));
    tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds             

    % define clipping ratios based on optimization results
    [tx.rclip, rx.rclip] = select_clipping_ratio(sim.type, ofdm.CS, tx.modulator.fc);
    
    %%
    berc = zeros(1, sim.Navg);
    bere = zeros(1, sim.Navg);
    Ptxm = zeros(1, sim.Navg);
    Ptxe = zeros(1, sim.Navg);
    for kk = 1:sim.Navg
        [ber, Ptxm(kk), Ptxe(kk)] = dc_ofdm(ofdm, tx, fiber, rx, sim);
   
        berc(kk) = ber.count;
        bere(kk) = ber.est;
    end
    
    % format some important variables
    PtxdBm(k) = 10*log10(mean(Ptxm)/1e-3);
    PtxedBm(k) = 10*log10(mean(Ptxe)/1e-3);
    bercount(k) = log10(mean(berc));
    berest(k) = log10(mean(bere));
end

%% Required power for ACO-OFDM AWGN (with ideal filters, without CP, without dc bias)
snr = fzero(@(x) berqam(ofdm.CS, x) - sim.BERtarget, 20);      % snr to achieve target BER;
snrl = 10^(snr/10);

K = 1 - 2*qfunc(tx.rclip);

Pnrx = snrl*ofdm.Ms*ofdm.Rs*rx.Sth/ofdm.Nc;
Pnawgn = K.^2*Pnrx./abs(rx.R*tx.kappa)^2*ones(1, ofdm.Nu/2);
Pawgn = tx.kappa*sqrt(2*sum(Pnawgn))*mean(tx.rclip);
PawgndBm = 10*log10(Pawgn/1e-3);

%% Required power for OOK AWGN
PookdBm = 10*log10(1/rx.R*qfuncinv(sim.BERtarget)*sqrt(rx.Sth*ofdm.Rb)/1e-3);

% Same as solving for square constellations
% qfuncinv((sim.BERtarget*log2(ofdm.CS)/4)/(1-1/sqrt(ofdm.CS)))/sqrt(3*pi*log2(ofdm.CS)/((ofdm.CS-1)*4*rx.Sth*ofdm.Rb))

% Power penalty with respect to NRZ-OOK assuming no bandwidth limitation.
% Oversampling penalty is included
PP_ofdm_dc = @ (M, r) 10*log10(r*sqrt(2*(M-1)./(3*log2(M)))); 

power_pen_ook_m = PtxdBm - PookdBm;
power_pen_ook_e = PtxedBm - PookdBm;

%% Figures
figure
subplot(211), hold on
plot(Fnl/1e9, bercount, '-sk', Fnl/1e9, berest, '-or')
plot(Fnl/1e9, log10(sim.BERtarget)*ones(size(Fnl)), 'k')
xlabel('Cut-off frequency (GHz)')
ylabel('log_{10}(BER)')
legend('BER counted', 'BER estimaded', 'Target BER')
axis([Fnl(1)/1e9 Fnl(end)/1e9 -4 -3])

subplot(212), hold on
plot(Fnl/1e9, power_pen_ook_m, '-sk', 'LineWidth', 2)
plot(Fnl/1e9, power_pen_ook_e, '-xb', 'LineWidth', 2)
plot(Fnl/1e9, (PawgndBm - PookdBm)*ones(size(Fnl)), '-r', 'LineWidth', 2)
xlabel('Cut-off frequency (GHz)', 'FontSize', 12)
ylabel('Optical Power Penalty (dB)', 'FontSize', 12)
legend('Measured', 'Estimated', 'AWGN')
