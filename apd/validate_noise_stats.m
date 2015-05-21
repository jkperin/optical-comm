%% clear, clc, close all
clear, clc, close all

format compact
addpath f

% Simulation
sim.BERtarget = 1.8e-4; % Target BER and SER
sim.Nsymb = 2^12;
sim.Mct = 1;    % Oversampling ratio for continuous-time simulation
sim.rin = 'off';
sim.levelSpacing = 'uniform';
sim.N0 = (20e-12)^2;

% M-PAM
mpam.M = 8;
mpam.Rb = 100e9;
mpam.bw = mpam.Rb/log2(mpam.M);
mpam.Rs = mpam.Rb/log2(mpam.M);

% TX
PrecdBm = -30;
Prec = 1e-3*10.^(PrecdBm/10)
tx.RIN = -150;
tx.lamb = 1310e-9;

% RX
% apd
% calss apd(GainBW, ka, Gain, noise_stats)
apd_gauss = apd(300e9, 0.09, 300e9/mpam.Rs, 'Gaussian');
apd_doubl = apd(300e9, 0.09, 300e9/mpam.Rs, 'DoublyStoch');

Pin = Prec*ones(1, sim.Nsymb);

out_gauss = apd_gauss.detect(Pin, 1/mpam.Rs);
out_doubl = apd_doubl.detect(Pin, 1/mpam.Rs);


figure, hist(out_gauss, 50);
figure, hist(out_doubl, 50);

[mean(out_gauss) mean(out_doubl)]
[var(out_gauss) var(out_doubl)]

pdf = apd_gauss.levels_pdf(unique(Pin), 1/mpam.Rs)


plot(pdf.I, pdf.p, pdf.I, pdf.p_gauss)
legend('accurate', 'gauss')

figure
plot(pdf.I, log10(pdf.p), pdf.I, log10(pdf.p_gauss))
legend('accurate', 'gauss')


