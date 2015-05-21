clear, clc, close all

format compact
addpath f

% Simulation
sim.BERtarget = 1.8e-4; % Target BER and SER
sim.Nsymb = 2^14;
sim.Mct = 1;    % Oversampling ratio for continuous-time simulation
sim.rin = 'off';
sim.levelSpacing = 'nonuniform';
sim.N0 = 0*(20e-12)^2;

% M-PAM
mpam.M = 8;
mpam.Rb = 100e9;
mpam.bw = mpam.Rb/log2(mpam.M);
mpam.Rs = mpam.Rb/log2(mpam.M);

% TX
PrecdBm = -30:2:-10;
Prec = 1e-3*10.^(PrecdBm/10);
tx.RIN = -150;
tx.lamb = 1310e-9;


% RX
% apd
% class apd(GainBW, ka, Gain, noise_stats)
% rx = apd(300e9, 0.09, 10, 'Gaussian');
rx = apd(300e9, 0.09, 10, 'DoublyStoch');

[a, b, a2, b2] = opt_level_spacing(0, mpam, rx, sim, sim.N0);

mpam1 = mpam;
mpam1.a = a;
mpam1.b = b;

mpam2 = mpam;
mpam2.a = a2;
mpam2.b = b2;


ber1 = zeros(size(Prec));
ber2 = zeros(size(Prec));
for k = 1:length(Prec)
    tx.Prec = Prec(k);
    
    ber1(k) = ber_mpam_imdd(mpam1, tx, rx, sim);
    
    ber2(k) = ber_mpam_imdd(mpam2, tx, rx, sim);
  
end


figure(1), hold on, box on, grid on
plot(PrecdBm, log10(ber1),  '-r', 'LineWidth', 1.5)
plot(PrecdBm, log10(ber2),  '-b', 'LineWidth', 1.5)
xlabel('P_{rec} (dBm)', 'FontSize', 12)
ylabel('log_{10}(BER)', 'FontSize', 12)
legend('Non-uniform spacing (numerical)', 'Non-uniform spacing (analytical)')
% title(sprintf('%d-PAM, %s', mpam.M, sim.levelSpacing))
axis([min(PrecdBm) -10 -7 0])



