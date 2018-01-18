%% Validate capacity in nonlinear regime
clear, close all

addpath ../data/
addpath ../
addpath ../f/
addpath ../../f/

S = load('../results/capacity_vs_pump_power_PdBm/capacity_vs_pump_power_EDF=principles_type3_pump=100mW_980nm_L=286_x_50km.mat');

E = S.Eopt{S.kopt};
Signal = S.Sopt{S.kopt};
Pump = S.Pump;
problem = S.problem;
problem.step_approx = @(x) 0.5*(tanh(2*x) + 1); % Smoothing factor = 2
problem.excess_noise = 1.5; % 1.2 for 980nm, 1.6 for 1480nm
problem.SwarmSize = min(100, 20*(Signal.N+1));

L0 = load('../../f/GN_model_coeff_spanLengthkm=50_l=0.mat');
Ln = load('../../f/GN_model_coeff_spanLengthkm=50_l=-1.mat');
Lp = load('../../f/GN_model_coeff_spanLengthkm=50_l=1.mat');

problem.nonlinear_coeff = {Ln.D, L0.D, Lp.D};

[lin_num, lin_approx] = capacity_linear_regime(E, Pump, Signal, problem)
[~, SE_lin_relx] = capacity_linear_regime_relaxed([E.L Signal.PdBm], E, Pump, Signal, problem);

[nlin_num, nllin_approx] = capacity_nonlinear_regime(E, Pump, Signal, problem)
% [~, SE_lin_relx] = capacity_linear_regime_relaxed([E.L Signal.PdBm], E, Pump, Signal, problem);

figure, hold on
plot(Signal.wavelength, lin_num.SE)
plot(Signal.wavelength, SE_lin_relx)
plot(Signal.wavelength, nlin_num.SE)