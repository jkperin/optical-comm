%% Validate class EDF and channels
clear, clc, close all

addpath ../
addpath ../f/

E = EDF(8, 'giles_ge:silicate');
% E = EDF(8, 'principles_type1');
E.plot('all');

spanAttdB = 9;
Pon = 5e-5;
Pump = Channels(980e-9, 50e-3, 'forward');
Signal = Channels(linspace(1530, 1565, 88)*1e-9, Pon, 'forward');
ASEf = Channels(Signal.wavelength, 0, 'forward');
ASEb = Channels(Signal.wavelength, 0, 'backward');
% E.plot('coefficients', Signal.wavelength)

Lopt = E.optimal_length(Pump, Signal, spanAttdB)

dlamb = Signal.wavelength(2)-Signal.wavelength(1);
df = E.c/Signal.wavelength(1)-E.c/Signal.wavelength(2);

GaindB_semi_analytical = E.semi_analytical_gain(Pump, Signal)
Pase_analytical = E.analytical_ASE_PSD(Pump, Signal)*df; % ASE power

[GaindB, Ppump_out, Psignal_out, Pase, sol] = E.two_level_system(Pump, Signal, ASEf, ASEb, df, 40)
% [n2, z] = E.metastable_level_population(sol, Signal, Pump, ASEf, true);


figure, plot(Signal.wavelength*1e9, GaindB, Signal.wavelength*1e9, GaindB_semi_analytical, '--')
figure, plot(Signal.wavelength*1e9, 10*log10(Pase/1e-3),...
    Signal.wavelength*1e9, 10*log10(Pase_analytical/1e-3), '--')
