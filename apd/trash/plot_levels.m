clear, clc, close all

addpath f

% Simulation
sim.BERtarget = 1.8e-4; % Target BER and SER
sim.Nsymb = 2^18;
sim.shot = 'on';
sim.rin = 'on';
sim.levelSpacing = 'uniform';

% M-PAM
mpam.M = 16;
mpam.Rb = 100e9;
mpam.bw = mpam.Rb/log2(mpam.M);

% TX
PrecdBm = 0;
Prec = 1e-3*10.^(PrecdBm/10);
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


[a, b] = calcOptLevelSpacing(mpam, tx, apd, sim);

levels = 1:16;

figure, hold on, box on, grid on
stem(levels, a/a(2), 'LineWidth', 1.5);
stem(levels, levels, 'r', 'LineWidth', 1.5)
xlabel('Level', 'FontSize', 12)
ylabel('Normalized Amplitude', 'FontSize', 12)
legend('Non-Uniform Spacing', 'Uniform Spacing', 'Location', 'NorthWest')
axis([levels(1) levels(end) 0 25])
set(gca, 'xtick', 1:16)
set(gca, 'ytick', 0:2.5:25)
