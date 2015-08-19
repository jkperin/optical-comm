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
tx.PtxdBm = -18:-4;
   
tx.lamb = 1310e-9; % wavelength
tx.alpha = 0; % chirp parameter
tx.RIN = -140;  % dB/Hz
tx.rexdB = -5;  % extinction ratio in dB. Defined as Pmin/Pmax

% Modulator frequency response
tx.modulator.fc = 50e9; % modulator cut off frequency
tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds
tx.modulator.H = @(f) (1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2));
tx.modulator.h = @(t) [0*t(t < 0) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0))];

%% Fiber
fiber = fiber(); % fiber(L, att(lamb), D(lamb))

%% Receiver
rx.N0 = (30e-12).^2; % thermal noise psd
rx.Id = 10e-9; % dark current
rx.R = 1; % responsivity
pin = apd(0, 0, Inf, rx.R, rx.Id);
% Electric Lowpass Filter
% rx.elefilt = design_filter('bessel', 5, mpam.Rs/(sim.fs/2));
rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);
rx.antialiasing = design_filter('bessel', 4, mpam.Rs/(sim.fs/2));

rx.eq.type = 'Fixed TD-SR-LE';
rx.eq.ros = 2;
rx.eq.Ntaps = 15;
rx.eq.Ntrain = 10e3;
rx.eq.mu = 1e-2;

% Transmitted power
Ptx = dBm2Watt(tx.PtxdBm);

%% Design equalizer
eq = rx.eq;
eq.Ntaps = 9;

Grx = rx.elefilt;

Hch = tx.modulator.H(sim.f).*exp(1j*2*pi*sim.f*tx.modulator.grpdelay)...
    .*fiber.Hfiber(sim.f, tx).*Grx.H(sim.f/sim.fs);

%% MMSE Time-domain symbol-rate equalizer
n = -floor(eq.Ntaps/2)*sim.Mct:sim.Mct*floor(eq.Ntaps/2);
v = 1/sqrt(sim.Mct)*mpam.pshape(n + 0*(sim.Mct-1)/2);
h = 1/sim.fs*tx.modulator.h(n/sim.fs + tx.modulator.grpdelay);
p = conv(h, v, 'same');
% p = p/max(abs(p));
% Calculate Grx impulse response and remove group delay
grx = impz(Grx.num, Grx.den, n);
grx = interp1(n - Grx.grpdelay, grx, n, 'spline');

% x = conv(p, grx, 'same');
p = conv(p, grx, 'same');

% Downsample to rate ros x Rs
dn = sim.Mct/eq.ros;
nd = [n(1):dn:-dn 0:dn:n(end)];
pd = interp1(n, p, nd);
pd = pd/abs(sum(pd)); % normalize to unit gain at DC

pd = pd(find(nd == 0) + (-floor(eq.Ntaps/2):floor(eq.Ntaps/2)));
nd = nd(find(nd == 0) + (-floor(eq.Ntaps/2):floor(eq.Ntaps/2)));


% % xd = [0 0.5 1 0.5 0];
% 
% % Toeplitz matrix
% X = toeplitz([xd.'; zeros(eq.Ntaps-1, 1)], [xd(1) zeros(1, eq.Ntaps-1)]);
% 
% % Get correct taps from Toeplitz matrix
% X = X(ceil(size(X, 1)/2)+(-floor(eq.Ntaps/2):floor(eq.Ntaps/2)), 1:eq.ros:end);

% e = zeros((eq.Ntaps+1)/eq.ros, 1); 
% e((length(e)+1)/2) = 1;
% 
% for k = 1:length(Ptx)
%     tx.Ptx = Ptx(k);
%     
%     % Overall link gain
%     link_gain = pin.Gain*fiber.link_attenuation(tx.lamb)*pin.R;
% 
%     % Ajust levels to desired transmitted power and extinction ratio
%     mpam.adjust_levels(tx.Ptx, tx.rexdB);
%     Pmax = mpam.a(end); % used in the automatic gain control stage
% 
%     % Modulated PAM signal
%     dataTX = randi([0 mpam.M-1], 1, sim.Nsymb); % Random sequence
%     xt = mpam.mod(dataTX, sim.Mct);
% %     xt(1:sim.Mct*sim.Ndiscard) = 0; % zero sim.Ndiscard first symbols
% %     xt(end-sim.Mct*sim.Ndiscard+1:end) = 0; % zero sim.Ndiscard last symbbols
% 
%     % Generate optical signal
%     [~, Pt] = optical_modulator(xt, tx, sim);
% 
%     %% Detect and add noises
%     yt = pin.detect(Pt, sim.fs, 'gaussian');
% 
%     % Automatic gain control
%     % Pmax = 2*tx.Ptx/(1 + 10^(-abs(tx.rexdB)/10)); % calculated from mpam.a
%     yt = yt/(Pmax*link_gain); % just refer power values back to transmitter
%     mpam.norm_levels;
% 
%     %% Equalization
%     if isfield(rx, 'eq')
%         rx.eq.TrainSeq = dataTX;
%     else % otherwise only filter using rx.elefilt
%         rx.eq.type = 'None';
%     end
% 
%     % NSR
%     rx.eq.NSR = (rx.N0*sim.fs/2)/(mean(abs(mpam.a).^2)*(Pmax*link_gain)^2);
%     % Note: make NSR = 0 for Zero forcing equalizer
%     [ydf, rx.eq] = equalize(rx.eq.type, yt, mpam, tx, fiber, rx, sim);
%     
%     %% Equalizer
% %     yt = real(ifft(fft(yt).*ifftshift(rx.antialiasing.H(sim.f/sim.fs))));
%     yt = real(ifft(fft(yt).*ifftshift(rx.elefilt.H(sim.f/sim.fs))));
% 
%     if mod(sim.Mct/eq.ros, 2) == 0
%         yk = yt(floor(sim.Mct/2)+1:sim.Mct/eq.ros:end);
%         tk = sim.t(floor(sim.Mct/2)+1:sim.Mct/eq.ros:end);
%     else % if sim.Mct is not multiple of ros, then interpolate
%         yk = interp1(1:length(yt), yt, (floor(sim.Mct/2)+1:sim.Mct/eq.ros:length(yt)).', 'spline');
%         tk = interp1(1:length(yt), sim.t, (floor(sim.Mct/2)+1:sim.Mct/eq.ros:length(yt)).');
% 
% %             yk = resample(ytaa, sim.ros, sim.Mct);
% %             tk = resample(sim.t, sim.ros, sim.Mct);            
%     end
%     
%     % NSR = noise signal ratio
%     if isfield(rx.eq, 'NSR') && false
%         W = X*(((X' + rx.eq.NSR*eye(size(X')))*X)\e);
%     else
%         W = X*((X'*X)\e);
%     end           
% 
%     W = W/abs(sum(W));
%     
%     
%     %% MMSE equalizer
%     yd2 = filter(W, 1, yk);
%     yd2 = circshift(yd2, [-(eq.Ntaps-1)/2+1 0]); % remove delay of transversal filter
%     
%     yd = yd2(1:eq.ros:end);
%     
%     figure, hold on
%     plot(t, yt)
%     plot(tk, yk, 'o')
%     plot(tk(1:eq.ros:end), yk(1:eq.ros:end), '*')
%     plot(tk(1:eq.ros:end), yd, '*')
%     
%     % Symbols to be discard in BER calculation
%     if length(yd) < sim.Nsymb
%         yd = [yd; zeros(sim.Nsymb-length(yd), 1)];
%     end
%     
%     Ndiscard = sim.Ndiscard*[1 1];
%     if isfield(rx.eq, 'Ntrain') 
%         Ndiscard(1) = Ndiscard(1) + rx.eq.Ntrain;
%     end
%     if isfield(rx.eq, 'Ntaps')
%         Ndiscard = Ndiscard + rx.eq.Ntaps;
%     end
%     ndiscard = [1:Ndiscard(1) sim.Nsymb-Ndiscard(2):sim.Nsymb];
%     yd(ndiscard) = []; 
%     dataTX(ndiscard) = [];
% 
%     % Demodulate
%     dataRX = mpam.demod(yd);
% 
%     % True BER
%     [~, ber_f(k)] = biterr(dataRX, dataTX);
% end
% rx.eq = rmfield(rx.eq, 'NSR');
ber = apd_ber(mpam, tx, fiber, pin, rx, sim);

figure, hold on
% plot(tx.PtxdBm, log10(ber_f), '--o')
plot(tx.PtxdBm, log10(ber.count), '--s')
plot(tx.PtxdBm, log10(ber.gauss), '-')
plot(tx.PtxdBm, log10(ber.awgn), '-')
legend('This', 'APD-BER', 'Estimated', 'AWGN')
axis([tx.PtxdBm([1 end]) -8 0])

% figure, hold on
% [Hw, w] = freqz(W, 1);
% plot(rx.eq.f, abs(rx.eq.H).^2)
% plot(w/(2*pi), abs(Hw).^2)
% % legend('Equaliz

