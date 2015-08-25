%% Generate waveforms for PAM and OFDM signals
clear, clc, close all

addpath ../f/
addpath ../apd/
addpath ../apd/f/
addpath ../mpam
addpath ../ofdm
addpath ../ofdm/f/

% Modulator frequency response
tx.modulator.fc = 50e9; % modulator cut off frequency
tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
tx.modulator.h = @(t) [0*t(t < 0) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0))];
tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

sim.Mct = 211;
sim.shot = false;

tx.Ptx = 1e-3*10^(0/10);
tx.rexdB = -10;

rx.N0 = (10e-12)^2;

pin = apd(0, 0, Inf);
fiber = fiber(0);

%% M-PAM
mpam = PAM(4, 100e9, 'equally-spaced', @(n) double(n >= 0 & n < sim.Mct));
mpam.adjust_levels(tx.Ptx, tx.rexdB);

%% Time and frequency
sim.Nsymb = 2^14;
sim.N = sim.Nsymb*sim.Mct;

rise_time = 10e-12;
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'

dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

sim.t = t;
sim.f = f;

% 
dataTX = randi([0 mpam.M-1], [sim.Nsymb 1]);

xpam = mpam.mod(dataTX, sim.Mct);

xpamf = filter(1, [1 -exp(-rise_time*sim.fs)], xpam);

[~, xpamf] = optical_modulator(xpam, tx, sim);

xpamf = pin.detect(xpamf, sim.fs, 'gaussian', rx.N0);

% eyediagram(circshift(xpamf(:), [-8 0]), 2*sim.Mct)

eyePAM = commscope.eyediagram(...
    'SamplingFrequency', sim.fs, ...
    'SamplesPerSymbol', sim.Mct, ...
    'OperationMode', 'Real Signal');

update(eyePAM, xpamf/max(abs(xpamf)))   
cmap = jet(64);
% cmap(1,:) = [1 1 1];
plot(eyePAM, cmap);


%% OFDM
sim.Nsymb = 2^10;
ofdm = ofdm(64, 52, 64, 100e9, 3, 3);

sim.type = 'preemphasis';
sim.BERtarget = 1e-4;
sim.quantiz = true;
sim.ENOB = 6;

tx.rclip = 4;
rx.rclip = 4;

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

clip_level = tx.rclip*sqrt(sum(2*ofdm.Pn.*abs(Gdac).^2));

% Signal std at transmitter and receiver
sigtx = sqrt(2*sum(ofdm.Pn));

% Check if power is at reasonable levels
assert(sigtx < 1e-3, 'Power is too high: %.2G (mW)\nCannot improve performance by increasing power.\n', 1e3*sqrt(2*sum(ofdm.Pn)))

%% Generate OFDM signal
xt = ofdm.generate_signal(tx, sim);

%% Clipping
xtc = xt;
xtc(xt < -clip_level) = -clip_level;
xtc(xt > clip_level) = clip_level;

%% Add dc bias (dc_bias determines whether it's added full dc bias or clipped dc bias)
xtc = xtc + clip_level;

%% Adjust signal to get to desired power
rex = 10^(-abs(tx.rexdB)/10);  % defined as Pmin/Pmax    

xmean  = mean(xtc);   % measured

%% Adjust transmitted power
if isfield(tx, 'Ptx') % if Ptx was provided, then scale signal to desired Ptx      
    % Scale and add additional dc bias due to finite extinction ratio
    xtc = xtc*tx.Ptx*(1 - 2*rex)/xmean + 2*tx.Ptx*rex; 
else % just add additional dc bias due to finite extinction ratio
    xtc = xtc + 2*xmean*rex; 
end

%% Apply frequency response of the laser
[~, xofdm] = optical_modulator(xtc, tx, sim);

xofdm = pin.detect(xofdm, sim.fs, 'gaussian', rx.N0);

% eyediagram(circshift(xofdm(:), [-8 0]), 2*sim.Mct)

figure
n = 1:(ofdm.Nc + ofdm.Npre_os);
n = repmat(n(:), floor(sim.N/(ofdm.Nc + ofdm.Npre_os)), 1);

%
% cloudPlot(n, xofdm, [], [], [150 150]);

% Create an eye diagram object
eyeObj = commscope.eyediagram(...
    'SamplingFrequency', sim.fs, ...
    'SamplesPerSymbol', ofdm.Nc + ofdm.Npre_os, ...
    'OperationMode', 'Real Signal');

% Update the eye diagram object with the transmitted signal
% update(eyeObj, xpam/max(abs(xpam)));

% update(eyeObj, xofdm/max(abs(xofdm)));
% 

% eyeObj.PlotType = '3D Color';

% figure
xofdmn = xofdm/max(abs(xofdm));
clip = (mean(xofdm) + std(xofdm)*2.3*[-1 1])/max(abs(xofdmn));
xofdmnc = xofdmn;
xofdmnc(xofdmn < clip(1)) = clip(1);
xofdmnc(xofdmn > clip(2)) = clip(2);

update(eyeObj, xofdmnc);
axis([0 (2*(ofdm.Nc + ofdm.Npre_os)-1)/sim.fs 0 1])
% hold on
% plot([0 (2*(ofdm.Nc + ofdm.Npre_os)-1)/sim.fs], clip(1)*[1 1], 'w')
% plot([0 (2*(ofdm.Nc + ofdm.Npre_os)-1)/sim.fs], clip(2)*[1 1], 'w')

cmap = jet(64);
% cmap(1,:) = [1 1 1];
plot(eyeObj, cmap);

% figure
% reset(eyeObj)
% update(eyeObj, xofdmnc);




% Manage the figures
% managescattereyefig(hFig, eyeObj, 'right');
