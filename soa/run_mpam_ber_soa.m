clear, clc, close all

format compact
addpath f

% Constants
h = 6.62606957e-34;
q = 1.60217657e-19;
c = 299792458;

% Simulation
sim.BERtarget = 1.8e-4; % Target BER and SER
sim.Nsymb = 2^14;
sim.Mct = 1;    % Oversampling ratio for continuous-time simulation
sim.soa = true;
sim.shot = 'off';
sim.rin = 'off';
sim.levelSpacing = 'uniform';

% M-PAM
mpam.M = 8;
mpam.Rb = 100e9;
mpam.Rs = mpam.Rb/log2(mpam.M);
mpam.bw = mpam.Rb/log2(mpam.M);

% TX
PrecdBm = -30:0.5:-10;
Prec = 1e-3*10.^(PrecdBm/10);
tx.RIN = -150;
tx.lamb = 1310e-9;

% RX
% apd
rx.N0 = (30e-12)^2;
rx.R = 1;

% SOA
soa.Fn = 9; % dB
soa.Seq = @(G) (G - 1)*10^(soa.Fn/10)/2*(h*c/tx.lamb); % one-sided PSD
soa.G = 10^(5/10);

% Generate unipolar PAM signal
if strcmp(sim.levelSpacing, 'uniform')
    mpam.a = (0:2:2*(mpam.M-1)).';
    mpam.b = (1:2:(2*(mpam.M-1)-1)).';
elseif strcmp(sim.levelSpacing, 'nonuniform')
    [mpam.a, mpam.b] = calcOptLevelSpacingAmpl(mpam, tx, soa, rx, sim);   
else
    error('Invalide Option!')
end

% SOA
soa.G = calcOptSOAGain(mpam, tx, soa, rx, sim, 100);

for k = 1:length(Prec)
    tx.Prec = Prec(k);
    
    mpam.b = mpam.b*tx.Prec/mean(mpam.a);
    mpam.a = mpam.a*tx.Prec/mean(mpam.a);
    
    % theoretical BER
    berm(k) = mpam_imdd_soa(mpam, tx, soa, rx, sim);
    
    bert(k) =  mpam_ber_soa(mpam, soa, rx);
    
%     bert2(k) = ber_mpam(mpam, tx, rx, sim);
    
end


figure(1), hold on, box on, grid on
plot(PrecdBm, log10(berm),  '-b', 'LineWidth', 1.5)
plot(PrecdBm, log10(bert),  ':r', 'LineWidth', 1.5)
% plot(PrecdBm, log10(bert2),  '--r', 'LineWidth', 1.5)
xlabel('P_{rec} (dBm)', 'FontSize', 12)
ylabel('log_{10}(BER)', 'FontSize', 12)
legend('Doubly-Stochastic', 'Gaussian')
title(sprintf('%d-PAM, %s', mpam.M, sim.levelSpacing))
axis([min(PrecdBm) max(PrecdBm) -7 0])



