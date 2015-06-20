%% Run transmission of dc-ofdm and check if with the selected clipping ratio 
%% is it possible to achieve target BER. Normally, when the error is higher
%% than the expected it's because the clippinig ratio is too small
% clear, clc, close all

format compact

addpath f           % functions path

rng('default')      % initiate default random number generator
rng('shuffle');     % Reinitialize the random number generator used by rand, randi, and randn with a seed based on the current time

%% Parameters to change
tx.channel = -3;                    % channels

L = (0:1:5)*1e3;                   % fiber length

tx.dlamb = 20e-9;                   % WDM channel spacing (m)

% sim.type = 'palloc';                % type of power allocation

% ofdm.CS = 16;                       % constellation size

% sim.ENOB = 5;                       % Effective number of bits for DAC and ADC (only used if sim.quantiz is true)

tx.fnl = 30e9;                      % modulator cut off frequency

tx.alpha = 2;                     % modulator chirp fac

sim.quantiz = true;                 % include quantization at both transmitter and receiver

sim.Nsymb = 2^14;                   % number of OFDM symbols

% sim.shot = true;                    % Include shot noise?
%     
% sim.RIN = true;                     % Include intensity noise?

tx.rexdB = 10;                      % exctinction ratio in dB

tx.RIN = -150;                      % RIN in dB/Hz. Only used if sim.RIN is true

verbose = false;                    % Show all plots? It'll slow donw simulation

% sim.Navg = 3;                     % Number of noise realizations

%% Simulation
% Additional parameters
imp_resolution = 1e-9;  % Minimum resolution of the impulse response of the filters
                        % the last value of the truncated, normalized
                        % impulse response must be smaller than this value.
                        

%%
%% Simulation parameters
sim.save = false;                   % save important data from simulation on file
% sim.type = 'preemphasis';           % type of power allocation
sim.Pb = 1.8e-4;                    % Target BER
% sim.Nsymb = 2^14;                   % number of OFDM symbols
sim.Mct = 8;                       % oversampling rate for emulation of continuous time 
                                    % Note: Mct has to be at least 4 so DT filter design approximates CT well enough./
                                                                      
%% Modulation parameters 
ofdm.ofdm = 'dc_ofdm';             % ofdm type
ofdm.full_dc = false;
ofdm.Nc = 64;                       % total number of subcarriers = FFT size (must be even)
ofdm.Nu = 52;                       % number of nonzero subcarriers (including complex conjugate)
% ofdm.CS = 16;                       % constellation size
ofdm.Rb = 107e9;                    % bit rate (total over two polarizations) (b/s)

ofdm.Ms = ofdm.Nc/ofdm.Nu;       % oversampling ratio
ofdm.Rs = 2*ofdm.Rb/log2(ofdm.CS);   % Symbol rate  
ofdm.B = ofdm.Nu/2*log2(ofdm.CS);    % Total number of bits per symbol

%% Fiber
fiber.att = 0.35;        % attenuation in dB/km
fiber.S0 = 0.092*1e3;    % dispersion slope (in s/m^3)
fiber.lamb0 = 1310e-9;   % zero-dispersion wavelength (also used to convert D to beta2)

tx.lamb = fiber.lamb0 + tx.channel*tx.dlamb;

fiber.D = @(lamb) fiber.S0/4*(lamb - fiber.lamb0^4./(lamb.^3)); % Dispersion curve

%% Transmitter parameters 
% tx.fnl = 35e9;                                  % laser cut off frequency
tx.Hl = @(f) 1./(1 + 2*1j*f/tx.fnl - (f/tx.fnl).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
tx.hl = @(t) (2*pi*tx.fnl)^2*t(t >= 0).*exp(-2*pi*tx.fnl*t(t >= 0));
tx.hl_delay = 2/(2*pi*tx.fnl);          % group delay of second-order filter in seconds

tx.kappa = 1;                         % current to optical power conversion (dc slope)

% tx.alpha = 2;                       % modulator chirp

% tx.dlamb = 20e-9;                      % WDM channel spacing (m)

% Transmitter filter (ZOH + filter)
tx.filter = design_filter('bessel', 5, 1/(ofdm.Ms*sim.Mct));

% Convolve with ZOH
bzoh = ones(1, sim.Mct);
tx.filter.num = conv(tx.filter.num, bzoh);
tx.filter.grpdelay = grpdelay(tx.filter.num, tx.filter.den, 1);
tx.filter.H = @(f) freqz(tx.filter.num, tx.filter.den, 2*pi*f).*exp(1j*2*pi*f*tx.filter.grpdelay);
tx.filter.noisebw = @(fs) noisebw(tx.filter.num, tx.filter.den, 2^15, fs); % equivalent two-sided noise bandwidth over larger number of samples (2^15) given a sampling frequency fs 

%% Receiver parameters
rx.R = 1;                           % responsivity
rx.NEP = 30e-12;                    % Noise equivalent power of the TIA at the receiver (A/sqrt(Hz))
rx.Sth = rx.R^2*rx.NEP^2/2;            % two-sided psd of thermal noise at the receiver (Sth = N0/2)

% Antialiasing filter
rx.filter = design_filter('gaussian', 4, 1/(ofdm.Ms*sim.Mct));

%% Calculate cyclic prefix
[ofdm.Npre_os, ofdm.Nneg_os, ofdm.Npos_os] = cyclic_prefix(ofdm, tx, rx, sim);

%% Time and frequency scales
ofdm.fs = ofdm.Rs*(ofdm.Nc + ofdm.Npre_os)/ofdm.Nu;               % sampling rate (Hz)
ofdm.fsct = sim.Mct*ofdm.fs;                                      % sampling frequency to emulate continuous time (Hz)
ofdm.fc = ofdm.fs/ofdm.Nc*(1:ofdm.Nu/2);                          % frequency at which subcarriers are located       

dt = 1/ofdm.fsct;                                                 % time increment in emulating continuous time (s)
Ntot = sim.Mct*(ofdm.Nc + ofdm.Npre_os)*sim.Nsymb;                % total number of points simulated in continuous time
df = ofdm.fsct/Ntot;                                              % frequency increment in continuous time (Hz)                    

sim.t = dt*(0:Ntot-1);                                            % continuous time scale
sim.f = -ofdm.fsct/2:df:ofdm.fsct/2-df;                           % frequency range in continuous time                    


%% define clipping ratios based on optimization results
% this script creates the variables sim.rcliptx, sim.rcliprx, RCLIPTX and
% RCLIPRX based on optimized values and on whether or not quantization is
% on
select_clipping_ratio

% Initiliaze variables
PtxdBm = zeros(size(L));
PtxedBm = zeros(size(L));
bercount = zeros(size(L));
berest = zeros(size(L));
PookdBm = zeros(size(L));

for k = 1:length(L)
    
    fiber.L = L(k);                       % channel number. 0 corresponds to zero-wavelength dispersion.
    
    for kk = 1:sim.Navg
        % generate, propagate, and detect
        [ber, P] = dc_ofdm(ofdm, tx, fiber, rx, sim, verbose);

        % format some important variables
        Ptx(kk) = P.Ptx;
        Ptxe(kk) = P.Ptx_est;
        berc(kk) = ber.count;
        bere(kk) = ber.est;
    end
    
    PtxdBm(k) = 10*log10(mean(Ptx)/1e-3);
    PtxedBm(k) = 10*log10(mean(Ptxe)/1e-3);
    bercount(k) = log10(mean(berc));
    berest(k) = log10(mean(bere));
    
    %% Required power for OOK AWGN
    PookdBm(k) = 10*log10(1/rx.R*qfuncinv(sim.Pb)*sqrt(rx.Sth*ofdm.Rb)/1e-3) + fiber.att*fiber.L/1e3;

    % Same as solving for square constellations
    % qfuncinv((sim.Pb*log2(ofdm.CS)/4)/(1-1/sqrt(ofdm.CS)))/sqrt(3*pi*log2(ofdm.CS)/((ofdm.CS-1)*4*rx.Sth*ofdm.Rb))
    
end

power_pen_ook_m = PtxdBm - PookdBm;
power_pen_ook_e = PtxedBm - PookdBm;


%% Figures
figure, hold on, grid on
plot(L/1e3, bercount, '-sk', L/1e3, berest, '-or')
plot(L/1e3, log10(sim.Pb)*ones(size(L)), 'k')
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

if sim.save
    filepath = ['results/' ofdm.ofdm '/'];
    filename = [sim.type '_CS=' num2str(ofdm.CS) '_ENOB=' num2str(sim.ENOB)];
    
    sim = rmfield(sim, {'f', 't'});
    
    save([filepath filename], 'ofdm', 'tx', 'rx', 'sim', 'fiber', 'ber', 'P', 'L')
    
    %saveas(gca, [filepath filename], 'fig');
end