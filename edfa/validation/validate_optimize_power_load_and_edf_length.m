%% Validate max_channels_on.m
clear, close all

addpath ../data/
addpath ../
addpath ../f/
addpath ../../f/

E = EDF(10, 'principles_type3');
% E.plot('all')

df = 50e9;
dlamb = df2dlamb(df);
lamb = 1522e-9:dlamb:1582e-9;
L = 14.3e6;
SMF = fiber(50e3, @(lamb) 0.18, @(lamb) 0);
Namp = round(L/SMF.L);

% Pon = 6e-4; % for 1W pump
% Pon = 1e-4; % for 100mW pump
Pon = 1e-4; % for < 100mW pump
Signal = Channels(lamb, Pon, 'forward');
Pump = Channels(980e-9, 65e-3, 'forward');

[~, spanAttdB] = SMF.link_attenuation(Signal.wavelength);

problem.Pon = Pon;
problem.spanAttdB = spanAttdB;
problem.df = df;
problem.Namp = Namp;
problem.step_approx = @(x) 0.5*(tanh(2*x) + 1); % Smoothing factor = 2
problem.excess_noise = 1.5; % 1.2 for 980nm, 1.6 for 1480nm
problem.SwarmSize = min(100, 20*(Signal.N+1));

% [Eopt_fmin, SignalOn_fmin] = optimize_power_load_and_edf_length('fminbnd', E, Pump, Signal, problem, true);
% Lopt2 = E.optimal_length(Pump, SignalOn_fmin, spanAttdB)
% GaindB = E.semi_analytical_gain(Pump, SignalOn_fmin);
% plot(Signal.wavelength*1e9, GaindB)
% drawnow

% [SEnum_fmin, SEapprox_fmin] = capacity_linear_regime(Eopt_fmin, Pump, SignalOn_fmin, spanAttdB, Namp, df)

% [Eopt_interp, SignalOn_interp] = optimize_power_load_and_edf_length('interp', E, Pump, Signal, problem, true);
% % Lopt2 = E.optimal_length(Pump, SignalOn_interp, spanAttdB)
% % GaindB = E.semi_analytical_gain(Pump, SignalOn_interp);
% % plot(Signal.wavelength*1e9, GaindB)
% % drawnow

[Eopt_pswarm, SignalOn_pswarm, exitflab, num, approx] = optimize_power_load_and_edf_length('particle swarm', E, Pump, Signal, problem, true);

% SE = capacity_linear_regime_relaxed([Eopt_pswarm.L SignalOn_pswarm.P], E, Pump, Signal, problem)
% [num, approx] = capacity_linear_regime(Eopt_pswarm, Pump, SignalOn_pswarm, problem)