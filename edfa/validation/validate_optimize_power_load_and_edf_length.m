%% Validate optimization of power load and EDF length using the particle swarm optimization
clear, clc, close all

addpath ../data/
addpath ../
addpath ../f/
addpath ../../f/

% E = EDF(10, 'giles_ge:silicate');
E = EDF(10, 'corning_type1');
% E = EDF(10, 'corning (NEW)');
% E.plot('all')

% Corning experimental data
% E = EDF(8, 'Corning (NEW)');
% E.core_radius = 1.19e-6;
% E.doping_radius = 0.98e-6; % Er3+ core radius. e.g., 1.2 um in [1, Table 1], 1.05um in [4, pg 156]
% E.rho0 = 6.68e18; % Er3+ concentraction (cm^3), e.g., 0.7e19 in [4, pg 156]
% E.NA = 0.23; % numerical aperture, e.g., 0.28 in [4, pg. 156]
% E.tau = 10e-3; % metastable lifetime in s        
% E.excess_loss = 0.28;

df = 50e9;  
dlamb = df2dlamb(df);
lamb = 1522e-9:dlamb:1582e-9;
% lamb = lamb(39:end); % > 1537 nm
% lamb = lamb(10:end); % > 15XX nm
L = 14.35e6;
SMF = fiber(50e3, @(l) 0.165*ones(size(l)), @(l) 20.4e-6*ones(size(l)));
SMF.gamma = 0.8e-3;
Namp = round(L/SMF.L);

% Pon = 6e-4; % for 1W pump
% Pon = 1e-4; % for 100mW pump
Pon =0.7e-4; % for < 100mW pump
Signal = Channels(lamb, Pon, 'forward');
Pump = Channels(980e-9, 30e-3, 'forward');

[~, spanAttdB] = SMF.link_attenuation(1550e-9); % same attenuation for all channels
spanAttdB = spanAttdB + 1.5;

problem.Pon = Pon;
problem.spanAttdB = spanAttdB;
problem.df = df;
problem.Gap = 10^(-1/10);
problem.Namp = Namp;
problem.step_approx = @(x) 0.5*(tanh(2*x) + 1); % Smoothing factor = 2
problem.diff_step_approx = @(x) sech(2*x).^2; % first derivative (used for computing gradient)
problem.excess_noise_correction = 1.6; %
problem.SwarmSize = min(100, 20*(Signal.N+1));
problem.nonlinearity = true;
S = load('../../f/GN_model_coeff_spanLengthkm=50km_Df=50GHz.mat');
problem.nonlinear_coeff = S.nonlinear_coeff;
problem.epsilon = 0.06; % From Fig. 17 of P. Poggiolini and I. Paper, “The GN Model
% of Non-Linear Propagation in Uncompensated Coherent Optical Systems,” 
% J. Lightw. Technol., vol. 30, no. 24, pp. 3857–3879, 2012.

options.AdaptationConstant = 0.1; 
options.FiniteDiffStepSize = 1e-6;
options.MaxIterations = 20;
options.AbsTol = 1e-3;
options.MinStep = 1e-4;
problem.saddle_free_newton.options = options;

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

[Eopt_pswarm, SignalOn_pswarm, exitflag_pswarm, num_pswarm, approx_pswarm] ... 
    = optimize_power_load_and_edf_length('particle swarm', E, Pump, Signal, problem, true);

% Performs a local gradient-based optimization after particle_swarm is done
% [Eopt_local, SignalOn_local, exitflag_local, num_local, approx_local] ... 
%     = optimize_power_load_and_edf_length('saddle-free newton', Eopt_pswarm, Pump, SignalOn_pswarm, problem, true);

