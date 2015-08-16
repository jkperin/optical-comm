%% Compare performance of M-PAM in 3 different scenarios
% 1. PIN receiver, no amplifier
% 2. SOA & PIN receiver
% 3. APD
clear, clc, close all

addpath ../f/
addpath ../apd/
addpath ../apd/f/

dBm2Watt = @(x) 1e-3*10.^(x/10);

%% Simulation parameters
sim.Nsymb = 2^14; % Number of symbols in montecarlo simulation
sim.Mct = 15;     % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
sim.BERtarget = 1e-4; 
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulations
sim.L = 2;

sim.RIN = false;
sim.shot = false;
sim.verbose = false;

%% M-PAM
mpam = PAM(4, 100e9, 'equally-spaced', @(n) double(n >= 0 & n < sim.Mct));

%% Time and frequency
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'

dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

sim.t = t;
sim.f = f;

%% Transmitter
tx.PtxdBm = -20:1:-4;
   
tx.lamb = 1310e-9; % wavelength
tx.alpha = 0; % chirp parameter
tx.RIN = -140;  % dB/Hz
tx.rexdB = -5;  % extinction ratio in dB. Defined as Pmin/Pmax

% Modulator frequency response
tx.modulator.fc = 30e9; % modulator cut off frequency
tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds
tx.modulator.H = @(f) (1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2));
tx.modulator.h = @(t) [0*t(t < 0) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0))];

%% Fiber
fiber = fiber(); % fiber(L, att(lamb), D(lamb))

%% Receiver
rx.N0 = (20e-12).^2; % thermal noise psd
rx.Id = 10e-9; % dark current
rx.R = 1; % responsivity
pin = apd(0, 0, Inf, rx.R, rx.Id);
% Electric Lowpass Filter
% rx.elefilt = design_filter('bessel', 5, mpam.Rs/(sim.fs/2));
rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);
rx.antialiasing = design_filter('bessel', 4, mpam.Rs/(sim.fs/2));

rx.eq.type = 'Fixed TD-SR-LE';
% rx.eq.ros = 2;
rx.eq.Ntaps = 11;
rx.eq.Ntrain = 2e3;
rx.eq.mu = 1e-2;

% Transmitted power
Ptx = dBm2Watt(tx.PtxdBm);

%% Design equalizer
n = -rx.eq.Ntaps*sim.Mct:sim.Mct*rx.eq.Ntaps;
v = mpam.pshape(n);
df = 1/201;
h = tx.modulator.h(n/sim.fs + tx.modulator.grpdelay);
h = h/max(h);
p = conv(h, v, 'same');
p = p/max(p);
x = conv(p, fliplr(p), 'same');
x = x/max(x);
plot(n, v, n, h, n, x)
nd = n(mod(n, sim.Mct) == 0);
xd = x(mod(n, sim.Mct) == 0);
hold on
stem(nd, xd)
X = toeplitz([xd.'; zeros(rx.eq.Ntaps-1, 1)], [xd(1) zeros(1, rx.eq.Ntaps-1)]);
e = zeros(rx.eq.Ntaps, 1);
e((rx.eq.Ntaps-1)/2) = 1;
cZF = X*((X'*X)\e);
cZF = cZF/abs(sum(cZF));

Gtx = design_filter('matched', mpam.pshape, 1/sim.Mct); 

Hmatched = conj(Gtx.H(sim.f/sim.fs).*tx.modulator.H(sim.f).*...
    exp(1j*2*pi*sim.f*tx.modulator.grpdelay).*fiber.Hfiber(sim.f, tx));


for k = 1:length(Ptx)
    tx.Ptx = Ptx(k);
    
    % Overall link gain
    link_gain = pin.Gain*fiber.link_attenuation(tx.lamb)*pin.R;

    % Ajust levels to desired transmitted power and extinction ratio
    mpam.adjust_levels(tx.Ptx, tx.rexdB);
    Pmax = mpam.a(end); % used in the automatic gain control stage

    % Modulated PAM signal
    dataTX = randi([0 mpam.M-1], 1, sim.Nsymb); % Random sequence
    xt = mpam.mod(dataTX, sim.Mct);
    xt(1:sim.Mct*sim.Ndiscard) = 0; % zero sim.Ndiscard first symbols
    xt(end-sim.Mct*sim.Ndiscard+1:end) = 0; % zero sim.Ndiscard last symbbols

    % Generate optical signal
    [~, Pt] = optical_modulator(xt, tx, sim);

    %% Detect and add noises
    yt = pin.detect(Pt, sim.fs, 'gaussian', rx.N0);

    % Automatic gain control
    % Pmax = 2*tx.Ptx/(1 + 10^(-abs(tx.rexdB)/10)); % calculated from mpam.a
    yt = yt/(Pmax*link_gain); % just refer power values back to transmitter
    mpam.norm_levels;

    %% Equalization
    if isfield(rx, 'eq')
        rx.eq.TrainSeq = dataTX;
    else % otherwise only filter using rx.elefilt
        rx.eq.type = 'None';
    end

    % Equalize
    [yd2, rx.eq] = equalize(rx.eq.type, yt, mpam, tx, fiber, rx, sim);
    %% Time-domain symbol-rate equalizer
    % matchedfilt = design_filter('matched', @(t) conv(mpam.pshape(t), 1/sim.fs*tx.modulator.h(t/sim.fs), 'full') , 1/sim.Mct); 
    yt = real(ifft(fft(yt).*ifftshift(Hmatched)));

    yk = yt(floor(sim.Mct/2)+1:sim.Mct:end);
    
    yd = filter(cZF, 1, yk);
    yd = circshift(yd, [-(length(cZF)-1)/2+1 0]);

    % Symbols to be discard in BER calculation
    Ndiscard = sim.Ndiscard*[1 1];
    if isfield(rx.eq, 'Ntrain')
        Ndiscard(1) = Ndiscard(1) + rx.eq.Ntrain;
    end
    if isfield(rx.eq, 'Ntaps')
        Ndiscard = Ndiscard + rx.eq.Ntaps;
    end
    ndiscard = [1:Ndiscard(1) sim.Nsymb-Ndiscard(2):sim.Nsymb];
    yd(ndiscard) = []; 
    dataTX(ndiscard) = [];

    % Demodulate
    dataRX = mpam.demod(yd);

    % True BER
    [~, ber_f(k)] = biterr(dataRX, dataTX);
end

ber = apd_ber(mpam, tx, fiber, pin, rx, sim);

figure, hold on
plot(tx.PtxdBm, log10(ber_f), '--o')
plot(tx.PtxdBm, log10(ber.count), '--s')
plot(tx.PtxdBm, log10(ber.gauss), '-')
legend('BER', 'Count', 'Estimated')
axis([tx.PtxdBm([1 end]) -8 0])

figure, hold on
C = freqz(cZF, 1, rx.eq.f, mpam.Rs);
plot(rx.eq.f, abs(rx.eq.H).^2)
plot(rx.eq.f, abs(C).^2)

