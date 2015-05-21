%% Sampling rate vs modulator cutoff frequency

clear, clc, close all

format compact

addpath f           % functions path
rng('default')      % initiate default random number generator

%% Parameters to swipe
CS = 16;

Fnl = (10:5:50)*1e9;

% Additional parameters
imp_resolution = 1e-9;  % Minimum resolution of the impulse response of the filters
                        % the last value of the truncated, normalized
                        % impulse response must be smaller than this value.

Nc = 64;
Nu = 52;
Ms = Nc/Nu;

                        
%% Simulation parameters
sim.Mct = 8;                        % oversampling rate for emulation of continuous time 
                                    % Note: Mct has to be at least 4 so DT filter design approximates CT well enough./
                                    
%% Transmitter parameters 
tx.kappa = 1;                         % current to optical power conversion (dc slope)

% Interpolation type. Could be one of the following
% {'ideal', 'butter', 'cheby1', 'ellipt', 'gaussian', 'bessel', 'fir'}
% Ideal is the ideal FIR interpolator design by interp
% The others refer to the type of filter that is used after ZOH
tx.filter = 'bessel';
tx.filter_order = 5;                % filter order
tx.filter_cutoff = 1/Ms;       % Cut-off frequency normalized by ofdm.fs
tx.imp_length = 300;                % length of the impulse response of the filter

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
rx.filter = 'butter';
rx.filter_order = 5;                            % filter order
rx.filter_cutoff = 1/Ms;                   % Cut-off frequency normalized by ofdm.fs
rx.imp_length = 400;                            % length of the impulse response of the filter

[bfilt, afilt] = design_filter(rx.filter, rx.filter_order, rx.filter_cutoff, sim.Mct);

rx.gadc = impz(bfilt, afilt, rx.imp_length).';  % impulse response
rx.gadc = rx.gadc/sum(rx.gadc);                 % normalize so frequency response at 0 Hz is 0 dB
rx.gadc_delay = grpdelay(bfilt, afilt, 1);

% Check if number of points used is high enough to attain desired
% resolution
assert(abs(rx.gadc(end)) < imp_resolution, 'rx.gadc length is not high enough');

%% DC-OFDM
dcofdm.ofdm = 'dc_ofdm';             % ofdm type
dcofdm.Nc = Nc;                       % total number of subcarriers = FFT size (must be even)
dcofdm.Nu = Nu;                     % number of nonzero subcarriers (including complex conjugate)
dcofdm.CS = CS;                       % Constellation size (effective CS in case of variable bit loading)    
dcofdm.Rb = 107e9;                    % bit rate (total over two polarizations) (b/s)

dcofdm.Ms = dcofdm.Nc/dcofdm.Nu;       % oversampling ratio
dcofdm.Rs = 2*dcofdm.Rb/log2(dcofdm.CS);   % Symbol rate  
dcofdm.B = dcofdm.Nu/2*log2(dcofdm.CS);    % Total number of bits per symbol

%% ACO-OFDM
acoofdm.ofdm = 'aco_ofdm';             % ofdm type
acoofdm.Nc = Nc;                       % total number of subcarriers = FFT size (must be even)
acoofdm.Nu = Nu/2;                     % number of nonzero subcarriers (including complex conjugate)
acoofdm.CS = CS;                       % Constellation size (effective CS in case of variable bit loading)    
acoofdm.Rb = 107e9;                    % bit rate (total over two polarizations) (b/s)

acoofdm.Ms = acoofdm.Nc/acoofdm.Nu;       % oversampling ratio
acoofdm.Rs = 2*acoofdm.Rb/log2(acoofdm.CS);   % Symbol rate  
acoofdm.B = acoofdm.Nu/2*log2(acoofdm.CS);    % Total number of bits per symbol

%% Iterate simulation varying the cut-off frequency of the 2nd-order filter of the transmitter
% Note: cut-off frequency changes the cyclic prefix length and thus chnages
% the sampling rate. 

% Initiliaze variables
dcofdm.Fs = zeros(size(Fnl));
dcofdm.Ncp = zeros(size(Fnl));

acoofdm.Fs = zeros(size(Fnl));
acoofdm.Ncp = zeros(size(Fnl));

for k = 1:length(Fnl)
    % Transmitter parameters
    tx.fnl = Fnl(k);
    tx.Hl = @(f) 1./(1 + 2*1j*f/tx.fnl - (f/tx.fnl).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
    tx.hl = @(t) (2*pi*tx.fnl)^2*t(t >= 0).*exp(-2*pi*tx.fnl*t(t >= 0));
    tx.hl_delay = 2/(2*pi*tx.fnl);          % group delay of second-order filter in seconds

    %% DC-OFDM
    % Calculate cyclic prefix
    [dcofdm.Npre_os, dcofdm.Nneg_os, dcofdm.Npos_os] = cyclic_prefix(dcofdm, tx, rx, sim);
    
    dcofdm.Ncp(k) = dcofdm.Npre_os;

    % Time and frequency scales
    dcofdm.Fs(k) = dcofdm.Rs*(dcofdm.Nc + dcofdm.Npre_os)/dcofdm.Nu;               % sampling rate (Hz)
                               
      
    %% ACO-OFDM
    % Calculate cyclic prefix
    [acoofdm.Npre_os, acoofdm.Nneg_os, acoofdm.Npos_os] = cyclic_prefix(acoofdm, tx, rx, sim);
    
    acoofdm.Ncp(k) = acoofdm.Npre_os;

    % Time and frequency scales
    acoofdm.Fs(k) = acoofdm.Rs*(acoofdm.Nc + acoofdm.Npre_os)/acoofdm.Nu;               % sampling rate (Hz)
    
end

figure1 = figure('Color',[1 1 1]);
hold on
plot(Fnl/1e9, dcofdm.Fs/1e9, '-ok', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
plot(Fnl/1e9, acoofdm.Fs/1e9, '-sk', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])

legend('DC-OFDM', 'ACO-OFDM')

% Create xlabel
xlabel('Frequency (GHz)','FontSize',14);

% Create ylabel
ylabel('Sampling rate (GHz)','FontSize',14);

set(gca, 'FontSize', 14)
% axis([14.99 50 4])
grid on
box on


figure1 = figure('Color',[1 1 1]);
hold on
plot(Fnl/1e9, dcofdm.Ncp, '-ok', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
plot(Fnl/1e9, acoofdm.Ncp, '-sk', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])

legend('DC-OFDM', 'ACO-OFDM')

% Create xlabel
xlabel('Frequency (GHz)','FontSize',14);

% Create ylabel
ylabel('Cyclic Prefix Length','FontSize',14);

set(gca, 'FontSize', 14)
% axis([14.99 50 4])
grid on
box on


