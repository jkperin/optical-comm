%% Plot noise distribution at different levels after APD detection modeled as doubly stochastic noise

clear, clc, close all

addpath f

% Simulation
sim.BERtarget = 1.8e-4; % Target BER and SER
sim.Nsymb = 2^10;
sim.shot = 'on';
sim.rin = 'off';
sim.levelSpacing = 'uniform';

% M-PAM
mpam.M = 8;
mpam.Rb = 100e9;
mpam.bw = mpam.Rb/log2(mpam.M);

% TX
tx.RIN = -150;

% RX
% apd
apd.N0 = (30e-12)^2;
apd.R = 1;
apd.ka = 0.09;
apd.Gbw = 300e9;
apd.Fa = @(ka, G) ka*G + (1 - ka)*(2 - 1/G);
% apd.Gapd = calcOptAPDGain(mpam, tx, apd, sim, apd.Gbw/mpam.bw) % Calculate optimal APD gain for target BER
apd.Gapd = calcOptAPDGain(mpam, tx, apd, sim, 1e3) % Calculate optimal APD gain for target BER


% pin
pin = apd;
pin.Gapd = 1;
pin.ka = 0;  

%
PrecdBm = -20;
tx.Prec = 1e-3*10.^(PrecdBm/10);

Plevels = 0:mpam.M-1;
Plevels = Plevels*tx.Prec/mean(Plevels);

Npoints = 2^12;
Nbins = 50;
dt = 1/mpam.bw;
for k = 2:length(Plevels)
    Ilevels = apd_doubly_stochastic(Plevels(k)*ones(1, Npoints), dt, tx, apd);
%     Ilevels = Ilevels + sqrt(apd.N0*mpam.bw)*randn(size(Ilevels));
    [Ihist, xhist] = hist(Ilevels, Nbins);
    Ihist = Ihist/sum(Ihist);
    figure(1), hold on
    plot(xhist, Ihist)
end