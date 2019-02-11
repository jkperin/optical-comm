clear, clc, close all

addpath ../
addpath ../f/
addpath ../../f/

S = load('../results/capacity_vs_pump_power_EDF=corning_type1_pump=980nm_L=286_x_50km/capacity_vs_pump_power_EDF=corning_type1_pump=150mW_980nm_L=286_x_50km.mat');

E = S.nlin.E;
Signal = S.nlin.S;
Pump = S.Pump;
problem = S.problem;
problem.Gap = 10^(-1/10);
problem.step_approx = @(x) 0.5*(tanh(2*x) + 1); % Smoothing factor = 2
problem.diff_step_approx = @(x) sech(2*x).^2; % First derivative of step_approx
problem.excess_noise = E.analytical_excess_noise(Pump, Signal);
problem.excess_noise = problem.excess_noise*problem.excess_noise_correction;

% optimize_power_load_and_edf_length('none', E, Pump, Signal, problem, true);
% 
% options.AdaptationConstant = 0.1; 
% options.FiniteDiffStepSize = 1e-6;
% options.MaxIterations = 10;
% options.AbsTol = 1e-3;
% options.MinStep = 1e-3;
% problem.saddle_free_newton.options = options;
% optimize_power_load_and_edf_length('saddle-free newton', E, Pump, Signal, problem, true);

Signal.P(Signal.P == 0) = eps;

[GaindB, Psignal_out, Ppump_out, dGaindB] = E.semi_analytical_gain(Pump, Signal);

[dNL, NL] = GN_model_noise_gradient(Signal.P, problem.nonlinear_coeff);
dNL = dNL*problem.Namp^(1+problem.epsilon); % Scale gradient to "Namp" spans

% Multiply by gain, since dNL is originally computed with respect to 
% the launch power, while the optimzation is done with the input power
% to the amplifier
dNL = dNL.*(10.^(problem.spanAttdB/10)); 

clims = [min(dNL(:)) max(dNL(:))];

figure, imagesc(dGaindB)
xlabel('P_k')
ylabel('\partial G_n/\partial P_k')
colorbar
figure, imagesc(dNL)
xlabel('P_k')
ylabel('\partial NL_n/\partial P_k')
colorbar
% figure, hold on
% plot(Signal.lnm, Signal.PdBm)
% % plot(S.nlin_local.S.lnm, S.nlin_local.S.PdBm)
% % plot(SignalOn_local.lnm, SignalOn_local.PdBm)
% plot(Snz.lnm, Snz.PdBm)
% ylim([-20 -10])