%% Plots
clear, clc, close all

addpath ../f/
addpath f/

f = linspace(0, 50, 100)*1e9;
wavelength = 1380e-9;
Fiber = fiber(0);

%% Fiber frequency response, chirp = 0
alpha = 0;
DL = [5 50 100 150];

figure, hold on, box on
for k = 1:length(DL)
    Fiber.L = 1e-3*DL(k)/Fiber.D(wavelength);
    plot(f/1e9, abs(Fiber.Himdd(f, wavelength, alpha, 'large signal')).^2, 'LineWidth', 2, 'DisplayName', sprintf('%d ps/nm', round(1e6*Fiber.D(wavelength)*Fiber.L/1e3)))
end
xlabel('Frequency (GHz)', 'FontSize', 12)
ylabel('|H_{IM-DD}(f)|^2', 'FontSize', 12)
leg = legend('-DynamicLegend');
set(leg, 'FontSize', 12);
set(gca, 'FontSize', 12)
m = matlab2tikz(gca);
m.write('Himdd-0chirp.tex');

%% Fiber frequency response, chirp = -1
alpha = -1;
DL = [5 50 100 150];
figure, hold on, box on
for k = 1:length(DL)
    Fiber.L = 1e-3*DL(k)/Fiber.D(wavelength);
    plot(f/1e9, abs(Fiber.Himdd(f, wavelength, alpha, 'large signal')).^2, 'LineWidth', 2, 'DisplayName', sprintf('%d ps/nm', round(1e6*Fiber.D(wavelength)*Fiber.L/1e3)))
end
xlabel('Frequency (GHz)', 'FontSize', 12)
ylabel('|H_{IM-DD}(f)|^2', 'FontSize', 12)
leg = legend('-DynamicLegend');
set(leg, 'FontSize', 12);
set(gca, 'FontSize', 12)
m = matlab2tikz(gca);
m.write('Himdd--1chirp.tex');

%% Power allocation with chirp = 0, D = 50
alpha = -1;
DL = 150;
Fiber.L = 1e-3*DL/Fiber.D(wavelength);
OFDM = ofdm(128, 100, 16, 100e9, 'palloc');
OFDM.set_cyclic_prefix(5, 5);
Gch = Fiber.Himdd(OFDM.fc, wavelength, alpha, 'large signal');
OFDM.power_allocation(Gch, ones(size(Gch)), 1.8e-4, true);


% Power allocation

