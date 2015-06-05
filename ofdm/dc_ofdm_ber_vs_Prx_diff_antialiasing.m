%% Run transmission of dc-ofdm and check if with the selected clipping ratio 
%% is it possible to achieve target BER. Normally, when the error is higher
%% than the expected it's because the clippinig ratio is too small
clear, close all

format compact

addpath f           % functions path

rng('shuffle');     % Reinitialize the random number generator used by rand, randi, and randn with a seed based on the current time

%% Parameters to change
sim.ENOB = 6;                       % Effective number of bits for DAC and ADC (only used if sim.quantiz is true)
sim.quantiz = true;                 % include quantization at both transmitter and receiver
sim.shot = true;                    % Include shot noise?
sim.RIN = true;                     % Include intensity noise?


verbose = false;                    % Show all plots? It'll slow donw simulation

sim.Navg = 1;                     % Number of noise realizations

%% Simulation
% Additional parameters
imp_resolution = 1e-9;  % Minimum resolution of the impulse response of the filters
                        % the last value of the truncated, normalized
                        % impulse response must be smaller than this value.
                        

%%
%% Simulation parameters
sim.save = false;                   % save important data from simulation on file
sim.type = 'preemphasis';           % type of power allocation
sim.Pb = 1.8e-4;                    % Target BER
sim.Nsymb = 2^12;                   % number of OFDM symbols
sim.Mct = 4;                       % oversampling rate for emulation of continuous time 
                                    % Note: Mct has to be at least 4 so DT filter design approximates CT well enough./
                                                                      
%% Modulation parameters 
ofdm.ofdm = 'dc_ofdm';             % ofdm type
ofdm.full_dc = false;
ofdm.Nc = 64;                       % total number of subcarriers = FFT size (must be even)
ofdm.Nu = 52;                       % number of nonzero subcarriers (including complex conjugate)
ofdm.CS = 16;                       % constellation size
ofdm.Rb = 112e9;                    % bit rate (total over two polarizations) (b/s)

ofdm.Ms = ofdm.Nc/ofdm.Nu;       % oversampling ratio
ofdm.Rs = 2*ofdm.Rb/log2(ofdm.CS);   % Symbol rate  
ofdm.B = ofdm.Nu/2*log2(ofdm.CS);    % Total number of bits per symbol

%% Fiber
fiber.att = 0;        % attenuation in dB/km
fiber.S0 = 0.092*1e3;    % dispersion slope (in s/m^3)
fiber.lamb0 = 1310e-9;   % zero-dispersion wavelength (also used to convert D to beta2)
fiber.L = 0;

% tx.dlamb = 20e-9;                      % WDM channel spacing (m)
% tx.channel = 0;                    % channels
% tx.lamb = fiber.lamb0 + tx.channel*tx.dlamb;
tx.lamb= 1296.2e-9;

fiber.D = @(lamb) fiber.S0/4*(lamb - fiber.lamb0^4./(lamb.^3)); % Dispersion curve

%% Transmitter parameters 
tx.fnl = 20e9;                      % modulator cut off frequency
tx.Hl = @(f) 1./(1 + 2*1j*f/tx.fnl - (f/tx.fnl).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
tx.hl = @(t) (2*pi*tx.fnl)^2*t(t >= 0).*exp(-2*pi*tx.fnl*t(t >= 0));
tx.hl_delay = 2/(2*pi*tx.fnl);          % group delay of second-order filter in seconds

tx.kappa = 1;                         % current to optical power conversion (dc slope)

tx.alpha = 1;                       % modulator chirp

tx.rexdB = 10;                      % exctinction ratio in dB
tx.RIN = -160;                      % RIN in dB/Hz. Only used if sim.RIN is true


% Interpolation type. Could be one of the following
% {'ideal', 'butter', 'cheby1', 'ellipt', 'gaussian', 'bessel', 'fir'}
% Ideal is the ideal FIR interpolator design by interp
% The others refer to the type of filter that is used after ZOH
tx.filter = 'bessel';
tx.filter_order = 5;                % filter order
tx.filter_cutoff = 1/ofdm.Ms;       % Cut-off frequency normalized by ofdm.fs
tx.imp_length = 512;                % length of the impulse response of the filter

[bfilt, afilt] = design_filter(tx.filter, tx.filter_order, tx.filter_cutoff, sim.Mct);
                 
switch tx.filter
    case 'ideal'
        tx.gdac = bfilt/sum(bfilt);       
        tx.gdac_delay = grpdelay(bfilt, afilt, 1);
    otherwise
        bzoh = ones(1, sim.Mct);    % Grp delay = (Mct-1)/2
        bdac = conv(bfilt, bzoh);   % numerator of DAC transfer function (Gzoh x Gfilt) 
        adac = afilt;               % denominator of DAC transfer function (Gzoh x Gfilt) ZOH is FIR
        tx.gdac = impz(bdac, adac, tx.imp_length).';    

        tx.gdac = tx.gdac/sum(tx.gdac);
        tx.gdac_delay = grpdelay(bdac, adac, 1);    % calculates the group delay by approximating the filter as an FIR filter
                                                    % whose impulse response is given by tx.gdac   

        % Check if number of points used is high enough to attain desired
        % resolution
        assert(abs(tx.gdac(end)) < imp_resolution, 'tx.gdac length is not high enough');
end

%% Receiver parameters
rx.R = 1;                           % responsivity
rx.NEP = 20e-12;                    % Noise equivalent power of the TIA at the receiver (A/sqrt(Hz))
rx.Sth = 2*rx.R^2*rx.NEP^2/2;            % two-sided psd of thermal noise at the receiver (Sth = N0/2)

% Antialiasing filter
rx.filter = 'gaussian';
rx.filter_order = 4;                            % filter order
rx.filter_cutoff = 1/ofdm.Ms;                   % Cut-off frequency normalized by ofdm.fs
rx.imp_length = 512;                            % length of the impulse response of the filter

[bfilt, afilt] = design_filter(rx.filter, rx.filter_order, rx.filter_cutoff, sim.Mct);

rx.gadc = impz(bfilt, afilt, rx.imp_length).';  % impulse response
rx.gadc = rx.gadc/sum(rx.gadc);                 % normalize so frequency response at 0 Hz is 0 dB
rx.gadc_delay = grpdelay(bfilt, afilt, 1);

% Check if number of points used is high enough to attain desired
% resolution
assert(abs(rx.gadc(end)) < imp_resolution, 'rx.gadc length is not high enough');

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
% sim.rcliptx = 3;
% sim.rcliprx = 3;

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
        ber1 = dc_ofdm(ofdm, tx, fiber, rx, sim1, verbose);
        
        ber2 = dc_ofdm(ofdm, tx, fiber, rx, sim2, verbose);

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
plot(PrxdBm, berest1, '-ob', PrxdBm, bercount1, '--xb')
plot(PrxdBm, berest2, '-or', PrxdBm, bercount2, '--xr')
plot(oepaper.PrxdBm, log10(oepaper.BER), '-sm')
plot(oepaper.PrxdBmexp, log10(oepaper.BERexp), '-ok')
xlabel('P_{rx} (dBm)')
ylabel('log_{10}(BER)')
legend('Preemphasis (Analysis)', 'Preemphasis (Montecarlo)', ...
    'Palloc (Analysis)', 'Palloc (Montecarlo)', 'oe-23-2-1176 Simulation B2B',...
    'oe-23-2-1176 Experiment B2B');
axis([-10 -2 -5 0])
% saveas(gca, sprintf('comparison_with_oepaper_%dkm', fiber.L/1e3), 'png')