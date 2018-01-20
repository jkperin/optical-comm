%% Validate class EDF and channels
clear, clc, close all

addpath ../data/
addpath ../
addpath ../f/
addpath ../../f/

% E = EDF(25, 'giles_ge:silicate');
E = EDF(15, 'corning_edf');
E.plot('all');

spanAttdB = 9;
Pon = 1e-5;
Pump = Channels(980e-9, 60e-3, 'forward');
Signal = Channels(linspace(1530, 1565, 80)*1e-9, Pon, 'forward');
ASEf = Channels(Signal.wavelength, 0, 'forward');
ASEb = Channels(Signal.wavelength, 0, 'backward');
% E.plot('coefficients', Signal.wavelength)

% Lopt = E.optimal_length(Pump, Signal, spanAttdB)

dlamb = Signal.wavelength(2)-Signal.wavelength(1);
df = E.c/Signal.wavelength(1)-E.c/Signal.wavelength(2);

GaindB_semi_analytical = E.semi_analytical_gain(Pump, Signal);
nsp_correction = 1.2; % factor to correct theoretical nsp
Pase_analytical = E.analytical_ASE_PSD(Pump, Signal, nsp_correction)*df; % ASE power

SignalOut = Signal;
[GaindB, Ppump_out, SignalOut.P, Pase, sol] = E.two_level_system(Pump, Signal, ASEf, ASEb, df, 100, true)

[n2, z] = E.metastable_level_population(sol, Signal, Pump, ASEf, true);

figure, hold on, box on
plot(Signal.wavelength*1e9, GaindB, 'DisplayName', 'Numerical')
plot(Signal.wavelength*1e9, GaindB_semi_analytical, 'DisplayName', 'Semi-analytical')
xlabel('Wavelength (nm)')
ylabel('Gain (dB)')

figure, hold on, box on
plot(Signal.wavelength*1e9, Watt2dBm(Pase), 'DisplayName', 'Numerical')
plot(Signal.wavelength*1e9, Watt2dBm(Pase_analytical), 'DisplayName', 'Semi-analytical')
xlabel('Wavelength (nm)')
ylabel('ASE power (dBm)')

