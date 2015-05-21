%% BER as a function of the clipping ratio at the transmitter only. 
%% Quantization is not included. If quantization is included, then clipping
%% ratio at the receiver must also be provided.
% Vary (palloc, preemphasis, ENOB, and CS (16,64)
clear, clc, close all

format compact

addpath f           % functions path
rng('default')      % initiate default random number generator

%% Parameters to swipe
Fnl = 20e9;

RCLIP = 3:0.25:5;

verbose = false;

Nave = 2; %Number of repetitions

% Additional parameters
imp_resolution = 1e-9;  % Minimum resolution of the impulse response of the filters
                        % the last value of the truncated, normalized
                        % impulse response must be smaller than this value.

%%
%% Simulation parameters
sim.save = false;                   % save important data from simulation on file
sim.type = 'palloc';                % Type of compensation at the transmitter {'preemphasis', 'palloc'}
sim.Pb = 1.8e-4;                    % Target BER
sim.Nsymb = 2^14;                   % number of OFDM symbols
sim.Mct = 8;                        % oversampling rate for emulation of continuous time 
                                    % Note: Mct has to be at least 4 so DT filter design approximates CT well enough./
                                    
% sim.rcliptx = 8;                      % Clipping ratio of DC-OFDM
                                    
sim.quantiz = false;                % include quantization at both transmitter and receiver
sim.ENOB = 5;                       % effective number of bits
% sim.rcliprx = 8;                      % Clipping ratio of DC-OFDM

%% Modulation parameters 
ofdm.ofdm = 'dc_ofdm';             % ofdm type
ofdm.Nc = 64;                       % total number of subcarriers = FFT size (must be even)
ofdm.Nu = 52;                     % number of nonzero subcarriers (including complex conjugate)
ofdm.CS = 64;                       % Constellation size (effective CS in case of variable bit loading)    
ofdm.Rb = 107e9;                    % bit rate (total over two polarizations) (b/s)

ofdm.Ms = ofdm.Nc/ofdm.Nu;       % oversampling ratio
ofdm.Rs = 2*ofdm.Rb/log2(ofdm.CS);   % Symbol rate  
ofdm.B = ofdm.Nu/2*log2(ofdm.CS);    % Total number of bits per symbol

%% Transmitter parameters 
tx.kappa = 1;                         % current to optical power conversion (dc slope)

% Transmitter parameters
tx.fnl = Fnl;
tx.Hl = @(f) 1./(1 + 2*1j*f/tx.fnl - (f/tx.fnl).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
tx.hl = @(t) (2*pi*tx.fnl)^2*t(t >= 0).*exp(-2*pi*tx.fnl*t(t >= 0));
tx.hl_delay = 2/(2*pi*tx.fnl);          % group delay of second-order filter in seconds

% Interpolation type. Could be one of the following
% {'ideal', 'butter', 'cheby1', 'ellipt', 'gaussian', 'bessel', 'fir'}
% Ideal is the ideal FIR interpolator design by interp
% The others refer to the type of filter that is used after ZOH
tx.filter = 'bessel';
tx.filter_order = 5;                % filter order
tx.filter_cutoff = 1/ofdm.Ms;       % Cut-off frequency normalized by ofdm.fs
tx.imp_length = 256;                % length of the impulse response of the filter

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
rx.NEP = 30e-12;                    % Noise equivalent power of the TIA at the receiver (A/sqrt(Hz))
rx.Sth = rx.R^2*rx.NEP^2/2;            % two-sided psd of thermal noise at the receiver (Sth = N0/2)

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

% Calculate cyclic prefix
[ofdm.Npre_os, ofdm.Nneg_os, ofdm.Npos_os] = cyclic_prefix(ofdm, tx, rx, sim);

% Time and frequency scales
ofdm.fs = ofdm.Rs*(ofdm.Nc + ofdm.Npre_os)/ofdm.Nu;               % sampling rate (Hz)
ofdm.fsct = sim.Mct*ofdm.fs;                                      % sampling frequency to emulate continuous time (Hz)
ofdm.fc = ofdm.fs/ofdm.Nc*(1:ofdm.Nu/2);                          % frequency at which subcarriers are located       

dt = 1/ofdm.fsct;                                                 % time increment in emulating continuous time (s)
Ntot = sim.Mct*(ofdm.Nc + ofdm.Npre_os)*sim.Nsymb;                % total number of points simulated in continuous time
df = ofdm.fsct/Ntot;                                              % frequency increment in continuous time (Hz)                    

sim.t = dt*(0:Ntot-1);                                            % continuous time scale
sim.f = -ofdm.fsct/2:df:ofdm.fsct/2-df;                           % frequency range in continuous time                    
    
%% Iterate simulation varying the cut-off frequency of the 2nd-order filter of the transmitter
% Note: cut-off frequency changes the cyclic prefix length and thus chnages
% the sampling rate. 

% Initiliaze variables
PtxdBm = zeros(size(RCLIP));
bercount = zeros(size(RCLIP));
berest = zeros(size(RCLIP));

rng_seed = rng;

for k = 1:length(RCLIP)
    sim.rcliptx = RCLIP(k);                    % Clipping ratio of ACO-OFDM
     
    PtxdBmv = zeros(1, Nave);
    berc = zeros(1, Nave);
    bere = zeros(1, Nave);                               
    
    %%
    rng(rng_seed);
    for kk = 1:Nave
        [ber, P] = dc_ofdm(ofdm, tx, rx, sim, verbose);
        
        PtxdBmv(kk) = 10*log10(P.Ptx/1e-3);
        berc(kk) = log10(ber.count);
        bere(kk) = log10(ber.est);
    end
   
    % format some important variables
    PtxdBm(k) = mean(PtxdBmv);
    bercount(k) = mean(berc);
    berest(k) = mean(bere); 
end

%% Required power for ACO-OFDM AWGN (with ideal filters, without CP, without dc bias)
snr = fzero(@(x) berqam(ofdm.CS, x) - sim.Pb, 20);      % snr to achieve target BER;
snrl = 10^(snr/10);

Pnrx = snrl*ofdm.Ms*ofdm.Rs*rx.Sth/ofdm.Nc;
Pnawgn = (1-qfunc(sim.rcliptx))^2*Pnrx./abs(rx.R*tx.kappa)^2*ones(1, ofdm.Nu/2);
Pawgn = tx.kappa*sqrt(2*sum(Pnawgn))*(sim.rcliptx*(1 - qfunc(sim.rcliptx)) + 1/sqrt(2*pi)*exp(-sim.rcliptx^2/2));
PawgndBm = 10*log10(Pawgn/1e-3);

%% Required power for OOK AWGN
PookdBm = 10*log10(1/rx.R*qfuncinv(sim.Pb)*sqrt(rx.Sth*ofdm.Rb)/1e-3);

% Same as solving for square constellations
% qfuncinv((sim.Pb*log2(ofdm.CS)/4)/(1-1/sqrt(ofdm.CS)))/sqrt(3*pi*log2(ofdm.CS)/((ofdm.CS-1)*4*rx.Sth*ofdm.Rb))

% Power penalty with respect to NRZ-OOK assuming no bandwidth limitation.
% Oversampling penalty is included
PP_ofdm_dc = @ (M, r) 10*log10(r*sqrt(2*(M-1)./(3*log2(M)))); 

power_pen_ook = PtxdBm - PookdBm;

%% Figures
figure
plot(RCLIP, bercount, '-sk', RCLIP, berest, '-or')
hold on
plot(RCLIP, log10(sim.Pb)*ones(size(RCLIP)), 'k')
% plot(RCLIP, berint.', 'xk', 'MarkerSize', 6)
xlabel('Clipping ratio', 'FontSize', 12)
ylabel('log_{10}(BER)', 'FontSize', 12)
legend('BER counted', 'BER estimaded', 'Target BER')

figure
plot(RCLIP, power_pen_ook, '-sk', 'LineWidth', 2)
hold on
plot(RCLIP, PP_ofdm_dc(ofdm.CS, sim.rcliptx)*ones(size(Fnl)), '-r', 'LineWidth', 2)
xlabel('Clipping ratio', 'FontSize', 12)
ylabel('Optical Power Penalty (dB)', 'FontSize', 12)
legend('Power penalty OFDM', 'Power penalty OFDM in AWGN')

