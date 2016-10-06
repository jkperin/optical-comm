clear, clc, close all

addpath ../../f % class fiber path

c = 299792458;  % speed of light
Fspacing = 100e9; % Hz
L = 40; % km
receiverSensitivity = -35; % dBm
maxD = 50e-6; 
att = 0.35; % dB/km
CDpenalty = 8;
Pmax = 9.8; % dBm

lamb0 = 1310e-9;
f0 = c/lamb0;
Dlamb = c/f0 - c/(f0 + Fspacing)

%% 100 GHz spacing
disp('> 100-GHz spacing')
Nch = 10;
Disp = 0;
while all(abs(Disp)*L <= maxD) 
    [lambchs, Disp] = optimal_wdm_channels(Dlamb, Nch);
    Nch = Nch + 1;
end
Nch  = Nch - 1

Plaunch = Pmax - 10*log10(Nch)
Margin = Plaunch - CDpenalty - att*L - receiverSensitivity

[lambchs, Disp] = optimal_wdm_channels(Dlamb, Nch);
Fiber = fiber(L*1e3);
lamb = linspace(lambchs(1), lambchs(end));
figure, hold on, box on
plot(lamb*1e9, 1e6*Fiber.D(lamb)*Fiber.L/1e3, 'k', 'LineWidth', 2)
stem(lambchs*1e9, 1e6*Disp*Fiber.L/1e3, 'k', 'LineWidth', 2, 'MarkerFaceColor', 'w')
xlabel('Wavelength (nm)', 'FontSize', 12)
ylabel('Dispersion (ps/nm)', 'FontSize', 12)
set(gca, 'FontSize', 12)

disp('> 200-GHz spacing')
Nch = 10;
Disp = 0;
while all(abs(Disp)*L <= maxD) 
    [lambchs, Disp] = optimal_wdm_channels(2*Dlamb, Nch);
    Nch = Nch + 1;
end
Nch  = Nch - 1

Plaunch = Pmax - 10*log10(Nch)
Margin = Plaunch - CDpenalty - att*L - receiverSensitivity

[lambchs, Disp] = optimal_wdm_channels(2*Dlamb, Nch);
Fiber = fiber(L*1e3);
lamb = linspace(lambchs(1), lambchs(end));
figure, hold on, box on
plot(lamb*1e9, 1e6*Fiber.D(lamb)*Fiber.L/1e3, 'k', 'LineWidth', 2)
stem(lambchs*1e9, 1e6*Disp*Fiber.L/1e3, 'k', 'LineWidth', 2, 'MarkerFaceColor', 'w')
xlabel('Wavelength (nm)', 'FontSize', 12)
ylabel('Dispersion (ps/nm)', 'FontSize', 12)
set(gca, 'FontSize', 12)



