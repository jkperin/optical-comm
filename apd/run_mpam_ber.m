clear, 

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
PrecdBm = -30:0;
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


% pin
pin = apd;
pin.Gapd = 1;
pin.ka = 0;  

ber = zeros(size(Prec));
for k = 1:length(Prec)
    tx.Prec = Prec(k);

    % measured BER
    berm(k) = mpam_imdd_unamplified(mpam, tx, apd, sim);

    % theoretical BER
    [mpam.a, mpam.b] = calcOptLevelSpacing(mpam, tx, apd, sim);
    bert(k) = ber_mpam(mpam, tx, apd, sim);
       
end


figure(1), hold on, box on, grid on
plot(PrecdBm, log10(berm),  '-b', 'LineWidth', 1.5)
plot(PrecdBm, log10(bert),  ':r', 'LineWidth', 1.5)
xlabel('P_{rec} (dBm)', 'FontSize', 12)
ylabel('log_{10}(BER)', 'FontSize', 12)
legend('Measured', 'Theory')
title(sprintf('%d-PAM, %s', mpam.M, sim.levelSpacing))
axis([min(PrecdBm) max(PrecdBm) -7 0])



