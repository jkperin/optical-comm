%% Calculate BER of system with SOA including crosstalk penalty
clear, clc, close all

addpath ../mpam
addpath ../f % general functions
addpath f

dBm2W = @(x) 1e-3*10.^(x/10);
h = 6.62606957e-34; % Planck
q = 1.60217657e-19; % electron charge
c = 299792458;      % speed of light

%% SOA Parameters
carrier_lifetime = 200e-12;
PsatdBm = 10;
G0 = 10^(10/10); % unsaturated SOA gain
Psat = dBm2W(PsatdBm);
Mch = 1; % Number of WDM channels

% Determine Gain x Input power curve
PfitdBm = linspace(-50, PsatdBm);
Pfit = dBm2W(PfitdBm);
Gfit = [G0 zeros(size(Pfit))];

for k = 2:length(Pfit)
    [Gfit(k), ~, exitflag] = fzero(@(G) Pfit(k) -  Psat/(abs(G)-1)*log(G0/abs(G)), Gfit(k-1));
       
    Gfit(k) = abs(Gfit(k));
    
    if exitflag ~= 1
        exitflag
    end
end

Gfit = Gfit(2:end);

%% Simulation parameters
sim.Nsymb = 2^15; % Number of symbols in montecarlo simulation
sim.Mct = 15;      % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
sim.Me = 16;       % Number of used eigenvalues
sim.L = 2; % de Bruijin sub-sequence length (ISI symbol length)
sim.BERtarget = 1e-4; 
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

sim.polarizer = ~true;
sim.shot = ~true; % include shot noise in montecarlo simulation 
sim.RIN = ~true; % include RIN noise in montecarlo simulation
sim.verbose = false; % show stuff

% M-PAM
% M, Rb, leve_spacing, pshape
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
tx.PtxdBm = -24:-14;

tx.lamb = 1310e-9;
tx.RIN = -150;  % dB/Hz
tx.rexdB = -5;  % extinction ratio in dB. Defined as Pmin/Pmax

%% fiber
fiber = fiber(); % back-to-back

%% Receiver
rx.N0 = (20e-12).^2; % thermal noise psd
rx.Id = 10e-9; % dark current
rx.R = 1; % responsivity
% Electric Lowpass Filter
% rx.elefilt = design_filter('bessel', 5, 0.5*mpam.Rs/(sim.fs/2));
rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);
% Optical Bandpass Filter
rx.optfilt = design_filter('fbg', 4, 200e9/(sim.fs/2));

%% SOA
soa = soa(20, 7, 1310e-9, 20); % soa(GaindB, NF, lambda, maxGaindB)

% KLSE Fourier Series Expansion 
[rx.U_fourier, rx.D_fourier, rx.Fmax_fourier] = klse_fourier(rx, sim, sim.Mct*(mpam.M^sim.L + 2*sim.L)); 

% BER without crosstalk
ber = soa_ber(mpam, tx, fiber, soa, rx, sim);

%% Crosstalk
% Transmitted power
Ptx = dBm2W(tx.PtxdBm);

% Auxiliary variables
f = sim.f/sim.fs;
Ns = floor(carrier_lifetime*mpam.Rs); % number of symbols to be averaged

for k = 1:length(Ptx)
    tx.Ptx = Ptx(k);
    
    % Ajust levels to desired transmitted power and extinction ratio
    mpam = mpam.adjust_levels(tx.Ptx, tx.rexdB);
    Pmax = mpam.a(end);

    % Modulated PAM signal
    dataTX = randi([0 mpam.M-1], Mch, sim.Nsymb); % Random sequence
    
    for kk = 1:Mch
        xt(kk, :) = mpam.mod(dataTX(kk, :), sim.Mct);
    end
        
    xt(:, 1:sim.Mct*sim.Ndiscard) = 0; % zero sim.Ndiscard first symbols
    xt(:, end-sim.Mct*sim.Ndiscard+1:end) = 0; % zero sim.Ndiscard last symbbols

    % Generate optical signal
    for kk = 1:Mch
       [Et(kk,:), Pt(kk,:)] = optical_modulator(xt(kk,:), tx, sim);
    end
       
    %% Amplifier with crosstalk
    Gsoa = zeros(length(Et), 1);
    Pin = zeros(length(Et), 1);
    for nn = Ns*sim.Mct:sim.N
        Pin(nn) = sum(mean(Pt(:, nn-Ns*sim.Mct+1:nn), 2));
        
        Gsoa(nn) = interp1(Pfit, Gfit, Pin(nn), 'spline');
    end
%     Gsoa(:) = interp1(Pfit, Gfit, Mch*tx.Ptx, 'spline');
    
    % assumming Gain >> 1
    Ssp = (Gsoa - 1)*10^(soa.Fn/10)/2*(h*c/soa.lamb); % Agrawal 6.1.15 3rd edition
    N0 = 2*Ssp; % one-sided baseband equivalent of Ssp

    NN = length(Gsoa);
    % N0 is the psd per polarization
    w_x = sqrt(1/2*N0*sim.fs/2).*(randn(NN, 1) + 1j*randn(NN, 1));
    w_y = sqrt(1/2*N0*sim.fs/2).*(randn(NN, 1) + 1j*randn(NN, 1));
    % Note: soa.N0 is one-sided baseband equivalent of ASE PSD.
    % Thus we multiply by sim.fs/2

    et = [Et(1,:).'.*sqrt(Gsoa) + w_x, w_y]; 
           
    % Optical bandpass filter
    eo = [ifft(fft(et(:, 1)).*ifftshift(rx.optfilt.H(f))),...
        ifft(fft(et(:, 2)).*ifftshift(rx.optfilt.H(f)))];

    %% Direct detection and add thermal noise
    %% Shot noise
    if isfield(sim, 'shot') && sim.shot
        q = 1.60217657e-19;      % electron charge (C)

        % Instataneous received power considering only attenuation from the fiber   
        Sshot = 2*q*(rx.R*sum(abs(eo).^2,2)  + rx.Id);     % one-sided shot noise PSD

        % Frequency is divided by two because PSD is one-sided
        wshot = sqrt(Sshot*sim.fs/2).*randn(size(et));
    else 
        wshot = 0;
    end

    % Direct detection and add noises
    if isfield(sim, 'polarizer') && ~sim.polarizer 
        yt = abs(eo(:, 1)).^2 + abs(eo(:, 2)).^2;
    else % by default assumes that polarizer is used
        yt = abs(eo(:, 1)).^2;
    end
    yt = yt + wshot + sqrt(rx.N0*sim.fs/2)*randn(sim.N, 1);

    % Overall link gain
    link_gain = mean(Gsoa)*fiber.link_attenuation(tx.lamb)*rx.R;
    
    % Automatic gain control
    % Pmax = 2*tx.Ptx/(1 + 10^(-abs(tx.rexdB)/10)); % calculated from mpam.a
    yt = yt./(Pmax*link_gain); % just refer power values back to transmitter
    mpam = mpam.norm_levels;

    %% Equalization
    if isfield(rx, 'eq')
        rx.eq.TrainSeq = dataTX(1, :);
    else % otherwise only filter using rx.elefilt
        rx.eq.type = 'None';
    end

    % Equalize
    [yd, rx.eq] = equalize(rx.eq.type, yt, mpam, tx, fiber, rx, sim);

    % Symbols to be discard in BER calculation
    Ndiscard = max(sim.Ndiscard, Ns*sim.Mct)*[1 1];
    if isfield(rx.eq, 'Ntrain')
        Ndiscard(1) = Ndiscard(1) + rx.eq.Ntrain;
    end
    if isfield(rx.eq, 'Ntaps')
        Ndiscard = Ndiscard + rx.eq.Ntaps;
    end
    ndiscard = [1:Ndiscard(1) sim.Nsymb-Ndiscard(2):sim.Nsymb];
    yd(ndiscard) = []; 
    dataTX(:, ndiscard) = [];

    % Demodulate
    dataRX = mpam.demod(yd);

    % True BER
    [~, ber.crosstalk(k)] = biterr(dataRX, dataTX(1, :));
end

figure, hold on, grid on
plot(tx.PtxdBm, log10(ber.count), '-o')
plot(tx.PtxdBm, log10(ber.est)) % KLSE Fourier
% plot(tx.PtxdBm, log10(ber_klse_freq))
% plot(tx.PtxdBm, log10(ber.gauss))
plot(tx.PtxdBm, log10(ber.awgn), 'k')
plot(tx.PtxdBm, log10(ber.crosstalk), '-o')
legend('Monte Carlo', 'KLSE Fourier & Saddlepoint Approx',... %  'KLSE Frequency Domain & Saddlepoint Approx',...
        'AWGN approximation', 'Monte Carlo w/ Crosstalk', 'Location', 'SouthWest')
axis([tx.PtxdBm([1 end]) -8 0])
xlabel('Transmitted Power (dBm)')
ylabel('log_{10}(BER)')

% Gain x Input Power
figure, box on
plot(PfitdBm, 10*log10(Gfit))
xlabel('Input Power (dBm)')
ylabel('Amplifier Gain (dB)')

% Gain Saturation
figure, box on
plot(PfitdBm, 10*log10(Gfit) + PfitdBm, PfitdBm, PfitdBm + 10*log10(G0))
xlabel('Input Power (dBm)')
ylabel('Amplifier Gain (dB)')
