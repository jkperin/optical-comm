%% Power distribution with Mch channels averaged over Ns symbols

clear, close all

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
G0 = 1800; % small-power Gain in linear units

Psat = dBm2W(PsatdBm);

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

figure, box on
plot(PfitdBm, 10*log10(Gfit))
xlabel('Input Power (dBm)')
ylabel('Amplifier Gain (dB)')

%% 
% M-PAM
% M, Rb, leve_spacing, pshape
mpam = PAM(4, 100e9, 'equally-spaced', @(n) double(n >= 0 & n < sim.Mct));
mpam.adjust_levels(1e-3*10^(-20/10), -5);

Ns = floor(carrier_lifetime*mpam.Rs); % number of symbols to be averaged
Mch = 4; % Number of WDM channels

% 0 <= ki <= Mch x Ns 
% k1 + k2 + k3 + k4 = Mch x Ns
K = randi([0 Mch*Ns], [2^20 4]);

Ksum = K(:,1)+K(:,2)+K(:,3)+K(:,4);

K(Ksum ~= Mch*Ns, :) = [];

Ptot = (K(:, 1)*mpam.a(1) + K(:, 2)*mpam.a(2) + K(:, 3)*mpam.a(3) + K(:, 4)*mpam.a(4))/Ns;

figure, histfit(Ptot, 20, 'normal')

Gsoa = interp1(Pfit, Gfit, Ptot, 'spline');

figure, histfit(Gsoa, 20, 'normal')

% hist(Gsoa, 20)
% K = unique(K);

% Simulation parameters
sim.Nsymb = 2^10; % Number of symbols in montecarlo simulation
sim.Mct = 7;      % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
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
tx.PtxdBm = -24:1:-14;

tx.lamb = 1310e-9;
tx.RIN = -150;  % dB/Hz
tx.rexdB = -10;  % extinction ratio in dB. Defined as Pmin/Pmax

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

% ber = soa_ber(mpam, tx, fiber, soa, rx, sim);
  
%% Crosstalk
% Transmitted power
Ptx = dBm2W(tx.PtxdBm);
Mch = 4;
% Ns = 10;

for k = 1:length(Ptx)
    tx.Ptx = Ptx(k);
    
    % Normalized frequency
    f = sim.f/sim.fs;

    % Overall link gain
    link_gain = soa.Gain*fiber.link_attenuation(tx.lamb)*rx.R;

    % Ajust levels to desired transmitted power and extinction ratio
    mpam.adjust_levels(tx.Ptx, tx.rexdB);

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
       
    % Amplifier
    %% !! with crosstalk
    Gsoa = zeros(length(Et), 1);
    Pin = zeros(length(Et), 1);
    for nn = Ns*sim.Mct:sim.N
        Pin(nn) = sum(mean(Pt(:, nn-Ns*sim.Mct+1:nn), 2));
        
        Gsoa(nn) = interp1(Pfit, Gfit, Pin(nn), 'spline');
    end
    
    % assumming Gain >> 1
    Ssp = (Gsoa - 1)*10^(soa.Fn/10)/2*(h*c/soa.lamb); % Agrawal 6.1.15 3rd edition

    N0 = 0*2*Ssp; % one-sided baseband equivalent of Ssp

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
%     yt = yt + wshot + sqrt(rx.N0*sim.fs/2)*randn(sim.N, 1);

    % Electric low-pass filter
    yt = real(ifft(fft(yt).*ifftshift(rx.elefilt.H(f))));

    % Sample
    ix = (sim.Mct-1)/2+1:sim.Mct:length(yt); % sampling points
    yd = yt(ix);

    % Discard first and last sim.Ndiscard symbols
    ndiscard = [1:sim.Ndiscard+2*Ns sim.Nsymb-sim.Ndiscard-2*Ns+1:sim.Nsymb];
    yd(ndiscard) = []; 
    dataTX(:, ndiscard) = [];

    % Automatic gain control
    yd = yd/link_gain; % just refer power values back to transmitter

    % Demodulate
    dataRX = mpam.demod(yd);

    % True BER
    [~, ber.crosstalk(k)] = biterr(dataRX, dataTX(1, :));
end

figure(1), hold on, grid on
plot(tx.PtxdBm, log10(ber.count), '-o')
plot(tx.PtxdBm, log10(ber.est)) % KLSE Fourier
% plot(tx.PtxdBm, log10(ber_klse_freq))
plot(tx.PtxdBm, log10(ber.gauss))
plot(tx.PtxdBm, log10(ber.awgn))
plot(tx.PtxdBm, log10(ber.crosstalk))
legend('Monte Carlo', 'KLSE Fourier & Saddlepoint Approx',... %  'KLSE Frequency Domain & Saddlepoint Approx',...
        'Gaussian Approximation', 'AWGN approximation',...
        'With Crosstalk',...
    'Location', 'SouthWest')
axis([tx.PtxdBm([1 end]) -8 0])
xlabel('Transmitted Power (dBm)')
ylabel('log_{10}(BER)')
