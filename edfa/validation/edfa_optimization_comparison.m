%% Compare optimization to the results presented in "SDM for Power-efficient Undersea Transmission" by O. K. Sinkin, et al
clear, clc, close all

addpath ../
addpath ../f/
addpath ../../f/

E = EDF(10, 'principles_type3');

df = 33e9;
dlamb = df2dlamb(df);
lamb = 1530e-9:dlamb:1565e-9;
L = 312*46e3;
SMF = fiber(46e3, @(lamb) 0.211, @(lamb) 0); % Attenuation to obtain 9.7 dB per span as in reference
Namp = round(L/SMF.L);

Pon = dBm2Watt(12)/82; % EDFA output power was 12 dBm
Signal = Channels(lamb, Pon, 'forward');
Pump = Channels(980e-9, 60e-3, 'forward');

[~, spanAttdB] = SMF.link_attenuation(Signal.wavelength);

problem.Pon = Pon;
problem.spanAttdB = spanAttdB;
problem.df = df;
problem.Namp = Namp;

[Eopt_fmin, SignalOn_fmin] = optimize_power_load_and_edf_length('fminbnd', E, Pump, Signal, problem, true);
Lopt2 = E.optimal_length(Pump, SignalOn_fmin, spanAttdB)
GaindB = E.semi_analytical_gain(Pump, SignalOn_fmin);
plot(Signal.wavelength*1e9, GaindB)
drawnow

% [Eopt_interp, SignalOn_interp] = optimize_power_load_and_edf_length('interp', E, Pump, Signal, problem, true);
% Lopt2 = E.optimal_length(Pump, SignalOn_interp, spanAttdB)
% GaindB = E.semi_analytical_gain(Pump, SignalOn_interp);
% plot(Signal.wavelength*1e9, GaindB)
% drawnow

problem.Pon = 1e-4;

[Eopt_pswarm, SignalOn_pswarm] = optimize_power_load_and_edf_length('particle swarm', E, Pump, Signal, problem, true);
