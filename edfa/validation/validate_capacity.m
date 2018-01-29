%% Validate capacity functions in linear and nonlinear regime
clear, close all

addpath ../data/
addpath ../
addpath ../f/
addpath ../../f/

S = load('../results/capacity_vs_pump_power_PdBm/capacity_vs_pump_power_EDF=principles_type3_pump=60mW_980nm_L=286_x_50km.mat');

E = S.Eopt{S.kopt};
Signal = S.Sopt{S.kopt};
Pump = S.Pump;
problem = S.problem;
problem.step_approx = @(x) 0.5*(tanh(2*x) + 1); % Smoothing factor = 2
problem.diff_step_approx = @(x) sech(2*x).^2; % First derivative of step_approx 
problem.excess_noise = 1.5; % 1.2 for 980nm, 1.6 for 1480nm
problem.epsilon = 0.05;

S = load('../../f/GN_model_coeff_spanLengthkm=50.mat');
problem.nonlinear_coeff = S.nonlinear_coeff;

[lin_num, lin_approx] = capacity_linear_regime(E, Pump, Signal, problem);
[~, SE_lin_relx] = capacity_linear_regime_relaxed([E.L Signal.PdBm], E, Pump, Signal, problem);
[nlin_num, nllin_approx] = capacity_nonlinear_regime(E, Pump, Signal, problem);
[~, ~, SE_nlin_relx] = capacity_nonlinear_regime_relaxed([E.L Signal.PdBm], E, Pump, Signal, problem);

% Continue optimization with fmincon. Gradient test only passes if gain is
% considered constant during optimization
% options = optimoptions('fmincon', 'Display', 'iter', 'UseParallel', true,...
%     'CheckGradients', true, 'SpecifyObjectiveGradient', true, 'FiniteDifferenceType', 'central');
% 
% [X, val, exitflag] = fmincon(@(X) capacity_nonlinear_regime_relaxed([E.L X], E, Pump, Signal, problem), ...
%     Signal.PdBm, [], [], [], [], [0 -50*ones(1, Signal.N)], [30 -5+zeros(1, Signal.N)], [], options);

figure, hold on
plot(Signal.lnm, lin_num.SE, 'DisplayName', 'Linear')
plot(Signal.lnm, SE_lin_relx, '--', 'DisplayName', 'Linear relaxed')
plot(Signal.lnm, nlin_num.SE, 'DisplayName', 'Nonlinear')
plot(Signal.lnm, SE_nlin_relx, '--k', 'DisplayName', 'Nonlinear relaxed')
xlabel('Wavelength (nm)')
ylabel('Spectral efficiency (bits/s/Hz)')
legend('-dynamiclegend')

figure, hold on
plot(Signal.lnm, Watt2dBm(nlin_num.NL), 'DisplayName', 'Nonlinear noise')
plot(Signal.lnm, Watt2dBm(nlin_num.Pase), 'DisplayName', 'ASE')
xlabel('Wavelength (nm)')
ylabel('Noise (dBm)')
legend('-dynamiclegend')

% figure, hold on
% plot(Signal.lnm, Signal.PdBm, 'DisplayName', 'Optimized')
% plot(Signal.lnm, X, 'DisplayName', 'continued with fmincon')
% xlabel('Wavelength (nm)')
% ylabel('Power loading (dBm)')
% legend('-dynamiclegend')