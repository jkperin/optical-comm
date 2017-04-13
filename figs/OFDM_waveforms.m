%% Generate waveforms for PAM and OFDM signals
clear, clc, close all

% load Pout_ook

addpath ../f/
addpath ../apd/
addpath ../apd/f/
addpath ../mpam
addpath ../ofdm
addpath ../ofdm/f/

% Modulator frequency response
tx.modulator.fc = 70e9; % modulator cut off frequency
tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
tx.modulator.h = @(t) [0*t(t < 0) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0))];
tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

type = 'aco-ofdm'

sim.Mct = 15;
% dt1 = 1/16;
% t1 = 0:dt1:(length(Pout)-1)*dt1;
% dt2 = 1/sim.Mct;
% t2 = 0:dt2:t1(end);

sim.shot = false;

tx.Ptx = 8e-3;
tx.rexdB = -5;

rx.N0 = (10e-12)^2;

pin = apd(0, 0, Inf);
fiber = fiber(0);

%% OFDM
sim.Nsymb = 2^7;
ofdm = ofdm(64, 52, 64, 100e9, 3, 3);

sim.type = 'preemphasis';
sim.BERtarget = 1e-4;
sim.quantiz = true;
sim.ENOB = 6;

tx.rclip = 1.5;
rx.rclip = 1.5;

tx.lamb = 1310e-9;

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
rx.Sth = rx.N0/2;       % two-sided psd of thermal noise at the receiver (Sth = N0/2)

% Antialiasing filter
rx.filter = design_filter('gaussian', 4, 1/(ofdm.Ms*sim.Mct));       

%% Time and frequency scales
sim.fs = sim.Mct*ofdm.fs;                                     % sampling frequency to emulate continuous time (Hz)  

sim.N = sim.Mct*(ofdm.Nc + ofdm.Npre_os)*sim.Nsymb;           % total number of points simulated in continuous time
         
dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

sim.t = t;
sim.f = f;  

%% Power allocation
ofdm.power_allocation(tx, fiber, rx, sim);

% The group delay is removed from the frequency response of the ADC and DAC
Gdac = tx.filter.H(ofdm.fc/sim.fs);       

if strcmp(type, 'aco-ofdm')
    ofdm.Pn(mod(1:ofdm.Nu/2, 2) == 0) = 0; % set even subcarriers to zero
    
    clip_level = [0 tx.rclip*sqrt(sum(2*ofdm.Pn.*abs(Gdac).^2))];
else
    clip_level = tx.rclip*sqrt(sum(2*ofdm.Pn.*abs(Gdac).^2))*[1 1];
end

% Signal std at transmitter and receiver
sigtx = sqrt(2*sum(ofdm.Pn));

% Check if power is at reasonable levels
assert(sigtx < 1e-3, 'Power is too high: %.2G (mW)\nCannot improve performance by increasing power.\n', 1e3*sqrt(2*sum(ofdm.Pn)))

%% Generate OFDM signal
xt = ofdm.generate_signal(tx, sim);

%% Clipping
xtc = xt;
xtc(xt < -clip_level(1)) = -clip_level(1);
xtc(xt > clip_level(2)) = clip_level(2);

% Normalize and add DC bias
xt = xt/clip_level(2);
xtc = xtc/clip_level(2);
clip_level = clip_level/clip_level(2);

xt = xt + clip_level(1);
xtc = xtc + clip_level(1);

clip_level(2) = clip_level(2) + clip_level(1);
clip_level(1) = 0;

% %% Add dc bias (dc_bias determines whether it's added full dc bias or clipped dc bias)
% xtc = xtc + clip_level;

%% Adjust signal to get to desired power
% rex = 10^(-abs(tx.rexdB)/10);  % defined as Pmin/Pmax    
% 
% xmean  = mean(xtc);   % measured
% 
% %% Adjust transmitted power
% if isfield(tx, 'Ptx') % if Ptx was provided, then scale signal to desired Ptx      
%     % Scale and add additional dc bias due to finite extinction ratio
%     xtc = xtc*tx.Ptx*(1 - 2*rex)/xmean + 2*tx.Ptx*rex; 
% else % just add additional dc bias due to finite extinction ratio
%     xtc = xtc + 2*xmean*rex; 
% end

%% Apply frequency response of the laser
% [~, xofdm] = optical_modulator(xtc, tx, sim);

% xofdm = pin.detect(xofdm, sim.fs, 'gaussian', rx.N0);

% eyediagram(circshift(xofdm(:), [-8 0]), 2*sim.Mct)

yt = reshape(xt, (ofdm.Nc + ofdm.Npre_os)*sim.Mct, sim.Nsymb);      % signal + noise
ytc = reshape(xtc, (ofdm.Nc + ofdm.Npre_os)*sim.Mct, sim.Nsymb);      % signal + noise

for k = 3 % plot 5 symbols
    figure, hold on, box on
    n = (1:(ofdm.Nc + ofdm.Npre_os)*sim.Mct).';
    plot(n, yt(:, k), ':k')
    plot(n, ytc(:, k), '-k')
    plot(n([1 end]), -clip_level(1)*[1 1], 'b')
    plot(n([1 end]), clip_level(2)*[1 1], 'b')
end

filename = sprintf('%s-waveform.tex', type);
matlab2tikz(filename)

