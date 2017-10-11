%% Validate max_channels_on.m
clear

addpath ../
addpath ../f/
addpath ../../f/

E = EDF(10, 'principles_type3');

df = 50e9;
dlamb = df2dlamb(df);
lamb = 1530e-9:dlamb:1565e-9;
L = 14e6;
SMF = fiber(50e3, @(lamb) 0.18, @(lamb) 0);
Namp = round(L/SMF.L);

Pon = 2e-05;
Signal = Channels(lamb, Pon, 'forward');
Pump = Channels(1480e-9, 50e-3, 'forward');

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

% [SEnum_fmin, SEapprox_fmin] = capacity_linear_regime(Eopt_fmin, Pump, SignalOn_fmin, spanAttdB, Namp, df)

% [Eopt_interp, SignalOn_interp] = optimize_power_load_and_edf_length('interp', E, Pump, Signal, problem, true);
% % Lopt2 = E.optimal_length(Pump, SignalOn_interp, spanAttdB)
% % GaindB = E.semi_analytical_gain(Pump, SignalOn_interp);
% % plot(Signal.wavelength*1e9, GaindB)
% % drawnow

problem.Pon = 5e-4;

[Eopt_pswarm, SignalOn_pswarm] = optimize_power_load_and_edf_length('particle swarm', E, Pump, Signal, problem, true);
