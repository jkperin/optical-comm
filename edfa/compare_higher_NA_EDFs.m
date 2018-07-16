%% Compare high-NA EDF with regular Corning EDFs
clear, clc, close all

addpath data
addpath f/
addpath ../f/

% Common simulation parameters
L = 14.35e6; % overall system length
SMF = fiber(50e3, @(l) 0.165*ones(size(l)), @(l) 20.4e-6*ones(size(l)));
SMF.gamma = 0.8e-3;
Namp = round(L/SMF.L); % number of amplifiers
[~, spanAttdB] = SMF.link_attenuation(1550e-9); % same attenuation for all channels
spanAttdB = spanAttdB + 1.5; % add 1.5 dB of margin in each span

df = 50e9; % WDM channel spacing
dlamb = df2dlamb(df);
lamb = 1522e-9:dlamb:1582e-9;
lamb(lamb > 1570e-9) = []; % discard wavelengths above 1570 nm

% Initialize signal and pump
Pon =0.7e-4; % for < 100mW pump
Signal = Channels(lamb, Pon, 'forward');
Pump = Channels(980e-9, 60e-3, 'forward');

% Problem parameters
problem.Pon = Pon;
problem.spanAttdB = spanAttdB;
problem.df = df;
problem.Gap = 10^(-1/10);
problem.Namp = Namp;
problem.step_approx = @(x) 0.5*(tanh(2*x) + 1); % Smoothing factor = 2
problem.diff_step_approx = @(x) sech(2*x).^2; % first derivative (used for computing gradient)
problem.excess_noise_correction = 1.4; %
problem.SwarmSize = min(100, 20*(Signal.N+1));
problem.nonlinearity = true;
S = load('../f/GN_model_coeff_spanLengthkm=50km_Df=50GHz.mat');
problem.nonlinear_coeff = S.nonlinear_coeff;
problem.epsilon = 0.06; % From Fig. 17 of P. Poggiolini and I. Paper, “The GN Model
% of Non-Linear Propagation in Uncompensated Coherent Optical Systems,” 
% J. Lightw. Technol., vol. 30, no. 24, pp. 3857–3879, 2012.

% Options for Saddle-free newton optimization
options.AdaptationConstant = 0.1; 
options.FiniteDiffStepSize = 1e-6;
options.MaxIterations = 20;
options.AbsTol = 1e-3;
options.MinStep = 1e-4;
problem.saddle_free_newton.options = options;

%% Regular EDF
E = EDF(10, 'corning (NEW)');
E.excess_loss = 0.3;

[Eopt_pswarm, SignalOn_pswarm, exitflag_pswarm, num_pswarm, approx_pswarm] ... 
    = optimize_power_load_and_edf_length('particle swarm', E, Pump, Signal, problem, true);

% Performs a local gradient-based optimization after particle_swarm is done
% [Eopt_local, SignalOn_local, exitflag_local, num_local, approx_local] ... 
%     = optimize_power_load_and_edf_length('saddle-free newton', Eopt_pswarm, Pump, SignalOn_pswarm, problem, true);

%% Higher-NA fiber
E = EDF(10, 'corning high NA');
E.excess_loss = 0.3;

[highNA_Eopt_pswarm, highNA_SignalOn_pswarm, highNA_exitflag_pswarm, highNA_num_pswarm, highNA_approx_pswarm] ... 
    = optimize_power_load_and_edf_length('particle swarm', E, Pump, Signal, problem, true);

% Performs a local gradient-based optimization after particle_swarm is done
% [highNA_Eopt_local, highNA_SignalOn_local, highNA_exitflag_local, highNA_num_local, highNA_approx_local] ... 
%     = optimize_power_load_and_edf_length('saddle-free newton', Eopt_pswarm, Pump, SignalOn_pswarm, problem, true);

disp('- Optimal EDF length:')
fprintf('Corning (NEW) = %.2f m\nCorning High NA = %.2f m\n', Eopt_pswarm.L, highNA_Eopt_pswarm.L)

disp('- Total fiber capacity')
fprintf('Corning (NEW) = %.2f Tb/s\nCorning High NA = %.2f Tb/s\n', sum(num_pswarm.SE*df)/1e12, sum(highNA_num_pswarm.SE*df)/1e12)



