%% Validate class EDF and channels
clear, clc, close all

addpath ../data/
addpath ../
addpath ../f/
addpath ../../f/

% E = EDF(25, 'giles_ge:silicate');
E = EDF(15, 'principles_type3');
% E.plot('all');

spanAttdB = 9;
Pon = 1e-5;
Pump = Channels(980e-9, 100e-3, 'forward');
Signal = Channels(linspace(1530, 1565, 40)*1e-9, Pon, 'forward');
ASEf = Channels(Signal.wavelength, 0, 'forward');
ASEb = Channels(Signal.wavelength, 0, 'backward');
% E.plot('coefficients', Signal.wavelength)

% Lopt = E.optimal_length(Pump, Signal, spanAttdB)

dlamb = Signal.wavelength(2)-Signal.wavelength(1);
df = E.c/Signal.wavelength(1)-E.c/Signal.wavelength(2);

GaindB_semi_analytical = E.semi_analytical_gain(Pump, Signal);
Pase_analytical = E.analytical_ASE_PSD(Pump, Signal)*df; % ASE power

SignalOut = Signal;
[GaindB, Ppump_out, SignalOut.P, Pase, sol] = E.two_level_system(Pump, Signal, ASEf, ASEb, df, 100, true)

[n2, z] = E.metastable_level_population(sol, Signal, Pump, ASEf, true);

[~, a] = E.absorption_coeff(980e-9)
[PCE, PCEmax] = E.power_conversion_efficiency(Pump, Signal, SignalOut)

fprintf('Photon flux difference = %g\n', (sum(Signal.P./Signal.Ephoton) + Pump.P/Pump.Ephoton) - sum(SignalOut.P./SignalOut.Ephoton))
