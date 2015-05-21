clear, clc, close all

addpath f

% Target SER
sim.BERtarget = 1.8e-4;

% Electron charge
q = 1.60217657e-19;

% M-PAM
mpam.M = 4;
mpam.bw = 100e9/log2(mpam.M);

% TX
PrecdBm = -30:0;
Prec = 1e-3*10.^(PrecdBm/10);
tx.RIN = -150;

% RX
apd.N0 = (30e-12)^2;
apd.R = 1;
apd.Gapd = 10;
apd.ka = 0.09;
apd.Fa = apd.ka*apd.Gapd + (1 - apd.ka)*(2 - 1/apd.Gapd);

% simulation
sim.SERtarget = sim.BERtarget*log2(mpam.M);
xi = qfuncinv(sim.SERtarget*mpam.M/(2*(mpam.M-1)));

% Define parameters so that total noise variance can be written in the form
% sigma^2 = alpha + beta*P + gamma*P^2

[a, b] = calcOptLevelSpacing(mpam, tx, apd, sim);

a/a(2)

b/a(2)
    

[Gapdopt, Precopt, exitflag] = fminbnd(@(Gapd) mean(calcOptLevelSpacing(mpam, tx, apd, sim, Gapd)), 1, 20)

%% apd
apd.Gapd = Gapdopt;
apd.Fa = apd.ka*apd.Gapd + (1 - apd.ka)*(2 - 1/apd.Gapd);

[a, b] = calcOptLevelSpacing(mpam, tx, apd, sim);

Popt_apd = mean(a);

ser_apd = 2*(mpam.M-1)/mpam.M*qfunc(Prec/Popt_apd*xi);

%% pin
pin = apd;
pin.Gapd =1;
pin.ka = 0;
pin.Fa = 0;

[a, b] = calcOptLevelSpacing(mpam, tx, pin, sim);

Popt_pin = mean(a);

ser_pin = 2*(mpam.M-1)/mpam.M*qfunc(Prec/Popt_pin*xi);


figure, hold on, box on
plot(PrecdBm, log10(ser_apd), 'b')
plot(PrecdBm, log10(ser_pin), '--r')
xlabel('P_{rec} (dBm)', 'FontSize', 12)
ylabel('log_{10}(SER)', 'FontSize', 12)
legend('apd (Thermal + Shot + RIN)', 'pin (Thermal + Shot + RIN)')
axis([-30 0 -10 0])
title(sprintf('SER vs Prec for %d-PAM -- apd: Gain = %.2f, k_A = %.2f', mpam.M, apd.Gapd, apd.ka))


