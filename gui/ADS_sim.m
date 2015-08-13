clear, clc, close all

addpath data
addpath ..\f\
addpath ..\mpam\
addpath ..\apd\

dBm2Watt = @(x) 1e-3*10.^(x/10);

load SampleData4PAMlong

%% Matlab Simulation
% Simulation parameters
sim.Nsymb = Nsymb; % Number of symbols in montecarlo simulation
sim.Mct = Mct;     % Oversampling ratio to simulate continuous time (must be odd) 
sim.N = N; % number points in 'continuous-time' simulation

sim.BERtarget = 1e-4;

sim.shot = ~true; % include shot noise in montecarlo simulation (always included for pin and apd case)
sim.RIN = true; % include RIN noise in montecarlo simulation

sim.verbose = false; % show stuff
sim.Ndiscard = 8; % number of symbols to be discarded from the begning and end of the sequence

% Equalization
eq.ros = 2;
eq.Ntaps = 15;
eq.Ntrain = 5e3;
eq.mu = 1e-2;

% Time and frequency
sim.fs = mpam.Rs*Mct;  % sampling frequency in 'continuous-time'

dt = Ts;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df+0.1).';

sim.t = t;
sim.f = f;

% Transmitter
tx.PtxdBm = -4;
% tx.PtxdBm = 10*log10(8);
tx.lamb = 1310e-9;
tx.alpha = 0;
tx.RIN = -140;
tx.rexdB = -5;   % exctinction ratio in dB. Defined as Pmin/Pmax

% Modulator frequency response
tx.modulator.fc = 24e9;
tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
tx.modulator.h = @(t) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0));
tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

tx_ads = tx;
tx_ads = rmfield(tx_ads, 'modulator');

% Fiber
fiber = fiber(0, @(l) 0, @(l) 0);

% Receiver
rx.R = 1; % responsivity
rx.NEP = (20e-12);
rx.N0 = (rx.NEP).^2; % thermal noise psd
rx.Sth = rx.N0/2;
rx.Id = 10e-9; % dark current

pin = apd(0, 0, Inf, rx.R, rx.Id);

% Electric low-pass filter
rise_time = 5e-12;
tau = rise_time/2.197;
a = exp(-Ts/tau);
imp = impz(1-a, [1 -a]);

% Matched filter
% rx.elefilt = design_filter('matched', @(n) conv(mpam.pshape(n), imp), 1/sim.Mct);
Hmatched = design_filter('matched', mpam.pshape, 1/sim.Mct);
rx.elefilt = design_filter('bessel', 5, 0.7*mpam.Rs/(sim.fs/2));
% rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);

% AWGN results
% Auxiliary variables
Deltaf = Hmatched.noisebw(sim.fs)/2; % electric filter one-sided noise bandwidth
% function to calculate noise std
varTherm = rx.N0*Deltaf; % variance of thermal noise

if sim.RIN
    varRIN =  @(Plevel) 10^(tx.RIN/10)*Plevel.^2*Deltaf;
else
    varRIN = @(Plevel) 0;
end

noise_std = @(Plevel) sqrt(varTherm + varRIN(Plevel) + pin.var_shot(Plevel, Deltaf));

% Overall link gain
link_gain = fiber.link_attenuation(tx.lamb)*pin.R;

% Zeroed samples 
ndiscard = [1:(sim.Ndiscard+eq.Ntrain) (sim.Nsymb-sim.Ndiscard-eq.Ntrain):sim.Nsymb];

dataTXref = mpam.demod(xtref(sim.Mct/2+1:sim.Mct:length(xtref)));

% Equalization
eq.TrainSeq = dataTXref;

Ptx = dBm2Watt(tx.PtxdBm);
for k = 1:length(Ptx)
    tx.Ptx = Ptx(k);

	% Ajust levels to desired transmitted power and extinction ratio
    mpam.adjust_levels(tx.Ptx, tx.rexdB);
    % BER for AWGN channel
    ber.awgn(k) = mpam.ber_awgn(noise_std);
   
    xt = mpam.mod(dataTXref, sim.Mct);
    
%     xt = filter(1-a, [1 -a], xt);

    Ptads = Poutnf.'/mean(Poutnf)*tx.Ptx;
    
    xt([sim.Mct*sim.Ndiscard end-sim.Mct*sim.Ndiscard:end]) = 0; % zero sim.Ndiscard last symbbols
    Ptads([sim.Mct*sim.Ndiscard end-sim.Mct*sim.Ndiscard:end]) = 0; % zero sim.Ndiscard last symbbols

	% Generate optical signal
    [Et, ~] = optical_modulator(xt, tx, sim);
    [Etads, ~] = optical_modulator(Ptads, tx_ads, sim); % without frequency response, which is already included in ADS simulation

	% Fiber propagation
    [~, Pt] = fiber.linear_propagation(Et, sim.f, tx.lamb);
    [~, Ptads] = fiber.linear_propagation(Etads, sim.f, tx.lamb);

    % Detect and add noises
    yt = pin.detect(Pt, sim.fs, 'gaussian', rx.N0);
    ytads = pin.detect(Ptads, sim.fs, 'gaussian', rx.N0);
  
    % Automatic gain control
    Pmax = 2*tx.Ptx/(1 + 10^(-abs(tx.rexdB)/10));
    yt = yt/(Pmax*link_gain); % just refer power values back to transmitter
    ytads = ytads/(Pmax*link_gain);
    mpam.norm_levels;
    
    % Equalization
    eq_type = 'None';
    [yd, Heq] = equalize(eq_type, yt, mpam, tx, fiber, rx, eq, sim);
    ydads = equalize(eq_type, ytads, mpam, tx, fiber, rx, eq, sim);
   
    % Discard first and last sim.Ndiscard symbols
    yd(ndiscard) = []; 
    ydads(ndiscard) = [];
    dataTX = dataTXref;
    dataTX(ndiscard) = [];
    
    % Detect ADS signal
    if mpam.M == 2
        ydads = ydads - mean(ytads);
        dataRXads = zeros(size(ydads)).';
        dataRXads(ydads > 0) = 1;
    else        
        dataRXads = mpam.demod(ydads);
    end
    
    % Demodulate
    dataRX = mpam.demod(yd);
    
    % BER
    [~, ber.sim(k)] = biterr(dataRX, dataTX);
    [~, ber.ads(k)] = biterr(dataRXads, dataTX);
    ber.awgn2(k) = mpam.ber_awgn(@(P) sqrt(rx.N0/2)/(Pmax*link_gain)*sqrt(trapz(sim.f, abs(Heq.^2))));
    1;
end

% [H, w] = freqz(Heq(end:-1:1), 1);
% figure, hold on
% plot(2*mpam.Rs/1e9*w/(2*pi), abs(H).^2)

figure, hold on, grid on, box on
plot(tx.PtxdBm, log10(ber.sim), '-o')
plot(tx.PtxdBm, log10(ber.ads), '-o')
plot(tx.PtxdBm, log10(ber.awgn), '-')
plot(tx.PtxdBm, log10(ber.awgn2), '--')
xlabel('Transmitted Power (dBm)')
ylabel('log_{10}(BER)')
legend('Matlab Simulation', 'ADS-DFB model', 'AWGN', 'AWGN noise enhanced')
axis([tx.PtxdBm([1 end]) -5 0])
% 
figure
plot(sim.f/1e9, abs(Heq).^2)

% figure
% plot(sim.f/1e9, abs(Heq).^2)
% xlabel('Frequency (GHz)')
% ylabel('|H_{eq}(f)|^2')
% axis([0 2*mpam.Rs/1e9 0 1.2*max(abs(Heq).^2)])

% rise_time = 5e-12;
% tau = rise_time/2.197;
% a = exp(-Ts/tau);
% imp = impz(1-a, [1 -a]);
% 
% xtf = filter(1-a, [1 -a], xt);
% 
% figure, hold on
% plot(xt)
% plot(xtf)
% 
% eyediagram(xt, 2*Mct)
% eyediagram(xtf, 2*Mct)
