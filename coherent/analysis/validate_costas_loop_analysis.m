%% Validate small-signal analysis of Costas loop for DP-QPSK
clear, close all

addpath ../f/
addpath ../../f/

Qpsk = QAM(4, 2*112e9);
N = 2^16;
totalLineWidth = 2*200e3; % sum of linewidths from LO and transmitter laser
Delay = 500e-12; % group delay in s
Ncpr = 1; % number of polarizations used in CPR
R = 1; % Responsivity
Plo = dBm2Watt(15); % LO power
Prx = dBm2Watt(-30); % received power
p = 2; % number of polarizations used
Kd = Ncpr*R*sqrt(2*Plo*Prx)/p; 
sim.ModFormat = Qpsk;
sim.BERtarget = 1e-4;   
ros = 1;
fs = ros*Qpsk.Rs;
sim.fs = fs;
q = 1.60217662e-19;
delaySamples = max(round(Delay*sim.fs), 1)
Delay = delaySamples/sim.fs;
sim.Ncpr = Ncpr;

% Optimize EPLL parameters
csi = sqrt(2)/2;
% % wn = 2*pi*0.7e9; % 
% wn = optimizePLL(csi, Delay, totalLineWidth, Ncpr, sim, true);
PLL.nums = [2*csi*wn wn^2];
PLL.dens = [1 0 0]; % descending powers of s
[PLL.numz, PLL.denz] = impinvar(PLL.nums, PLL.dens, sim.fs);
fprintf('Loop filter fn: %.3f GHz\n', wn/(2*pi*1e9));

% Loop filter
LoopFilter = ClassFilter(PLL.numz, PLL.denz, sim.fs);

ReceiverFilter = ClassFilter('bessel', 5, 0.5*Qpsk.Rs/(sim.fs/2));

Laser = laser(1310e-9, Watt2dBm(Plo), -Inf, totalLineWidth, 0);


E = sqrt(Plo)*ones(1, N);
[E, phiPN] = Laser.addPhaseNosie(E, fs);

Ppd = Plo/(4*p);
varAWGN = (2*2*q*R*Ppd*fs/2)*2*Ncpr;
nAWGN = sqrt(varAWGN)*randn(size(phiPN));
nAWGN = ReceiverFilter.filter(nAWGN);
[var(nAWGN)/varAWGN, Qpsk.Rs/fs]

SNR = R*Prx/(p*q*Qpsk.Rs);

phiE = zeros(size(phiPN));
phiLO = zeros(size(phiPN));
phiIN = zeros(size(phiPN));
for t = delaySamples+1:N
    phiE(t) = phiPN(t) - phiLO(t-delaySamples);
    
    phiIN(t) = sin(phiE(t)) + nAWGN(t)/Kd;
    
    phiLO(t) = LoopFilter.filter(phiIN(t));
end

figure, hold on, box on
plot(phiPN)
plot(phiLO)

Nstart = 1e4;
Nend = N-Nstart;

w = 2*pi*linspace(-fs/2, fs/2, 1e3);
numFs = PLL.nums ;
denFs = [1 0]; % removes integrator due to oscillator
    
Fw = polyval(numFs, 1j*w)./polyval(denFs, 1j*w);
H1 = 1j*w + exp(-1j*w*Delay).*Fw;    
Sawgn = (2*q*R*Ppd)*2*Ncpr; % double-sided AWGN PSD
% varphiE = totalLineWidth*trapz(w, 1./abs(H1).^2) + Sawgn/(2*pi*Kd^2)*trapz(w, abs(Fw./H1).^2)
% varphiE = totalLineWidth*trapz(w, 1./abs(H1).^2) + 1/(2*pi*Ncpr*2*SNR*Qpsk.Rs)*trapz(w, abs(Fw./H1).^2)
[varphiE, nPN, nAWGN] = phase_error_variance(csi, wn, Ncpr, Delay, totalLineWidth, 10*log10(SNR), Qpsk.Rs, true)

disp('Phase error variance:')
fprintf('Theory: %g\nMeasured: %g\nTheory/Measured: %g\n', varphiE, var(phiE(Nstart:Nend)), varphiE/var(phiE(Nstart:Nend)))

fprintf('Contribution of phase noise vs AWGN on PLL phase error at receiver sensitivity (SNRdB = %.2f):\nPN/AWGN = %.3f\n', 10*log10(SNR),...
    totalLineWidth*trapz(w, 1./abs(H1).^2)/(1/(2*pi*Ncpr*2*SNR*Qpsk.Rs)*trapz(w, abs(Fw./H1).^2)));

% mismatch penalty
% Deltat = 100/Qpsk.Rs;
% 1/4*2*pi*totalLineWidth*Deltat*trapz(w, abs(w./H1).^2)/(2*pi*fs)

