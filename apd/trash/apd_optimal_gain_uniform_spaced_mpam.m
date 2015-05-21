clear, clc, close all

addpath f

% M-PAM
mpam.M = 8;
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

% Target BER and SER
sim.BERtarget = 1.8e-4;
sim.SERtarget = sim.BERtarget*log2(mpam.M);

% [Prec, ~, exitflag] = fzero(@(Prec) ser_uniform_mpam_tsr(mpam, tx, rx, rx.Gapd, 1e-3*10^(Prec/10)) - Ps, -20);

[Gapdopt, Precopt, exitflag] = fminbnd(@(Gapd) fzero(@(Prec) ser_uniform_mpam_tsr(mpam, tx, apd, Gapd, 1e-3*10^(Prec/10)) - sim.SERtarget, -20), 1, 20)

ser_uniform_mpam_tsr(mpam, tx, apd, Gapdopt, 1e-3*10^(Precopt/10))

for M = 1:20
    [Prec2(M), ~, exitflag] = fzero(@(Prec) ser_uniform_mpam_tsr(mpam, tx, apd, M, 1e-3*10^(Prec/10)) - sim.SERtarget, -20);
end

figure, hold on, box on
plot(1:20, Prec2)
plot(Gapdopt, Precopt, 'xr')
xlabel('apd Gain', 'FontSize', 12)
ylabel('P_{rec} (dBm)', 'FontSize', 12)
title(sprintf('P_{rec} x apd Gain to achieve BER = %.1E for %d-PAM (uniform spacing)', sim.BERtarget, mpam.M))

% SER for optimal gain
apd.Gapd = Gapdopt;

ser_opt = zeros(size(Prec));
for k = 1:length(Prec)
    tx.Prec = Prec(k);
    ser_opt(k) = ser_uniform_mpam_tsr(mpam, tx, apd);
end

Ps_opt = 2*(mpam.M-1)/mpam.M*qfunc(apd.R*apd.Gapd*Prec/(mpam.M-1)*1/sqrt(apd.N0*mpam.bw));

% SER for pin
pin = apd;
pin.Gapd = 1;
pin.ka = 0;

ser_pin = zeros(size(Prec));
for k = 1:length(Prec)
    tx.Prec = Prec(k);
    ser_pin(k) = ser_uniform_mpam_tsr(mpam, tx, pin);
end

Ps_pin = 2*(mpam.M-1)/mpam.M*qfunc(pin.R*pin.Gapd*Prec/(mpam.M-1)*1/sqrt(pin.N0*mpam.bw));

figure, hold on, box on
plot(PrecdBm, log10(ser_opt), 'b')
plot(PrecdBm, log10(Ps_opt), '--b')
plot(PrecdBm, log10(ser_pin), 'r')
plot(PrecdBm, log10(Ps_pin), '--r')
xlabel('P_{rec} (dBm)', 'FontSize', 12)
ylabel('log_{10}(SER)', 'FontSize', 12)
legend('apd (Thermal + Shot + RIN)', 'apd (Thermal)', 'pin (Thermal + Shot + RIN)', 'pin (Thermal)','Location', 'SouthWest')
axis([-30 0 -10 0])
title(sprintf('SER vs Prec for %d-PAM -- apd: Gain = %.2f, k_A = %.2f', mpam.M, apd.Gapd, apd.ka))



