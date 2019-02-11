%% Test gain gradient
function test_gain_gradient()
%% Validate gradient calculations
clear, clc, close all

addpath ../
addpath ../../f

S = load('../results/capacity_vs_pump_power_PdBm/capacity_vs_pump_power_EDF=principles_type3_pump=60mW_980nm_L=286_x_50km.mat');

E = S.Eopt{S.kopt};
Signal = S.Sopt{S.kopt};
Pump = S.Pump;
problem = S.problem;
problem.step_approx = @(x) 0.5*(tanh(2*x) + 1); % Smoothing factor = 2
problem.diff_step_approx = @(x) sech(2*x).^2; % First derivative of step_approx
problem.excess_noise = 1.5; % 1.2 for 980nm, 1.6 for 1480nm
problem.epsilon = 0.05;

options = optimoptions('fmincon', 'Display', 'iter', 'UseParallel', false,...
    'CheckGradients', true, 'SpecifyObjectiveGradient', true, 'FiniteDifferenceType', 'central',...
    'FiniteDifferenceStepSize', 1e-10);

Signal = Signal.sample(100:120);
P = Signal.P;
P(P == 0) = eps;
% [PdBm, val, exitflag] = fmincon(@(PdBm) objective(PdBm, E, Pump, Signal), ...
%     10*log10(P/1e-3), [], [], [], [], -20 + zeros(size(P)), zeros(size(P)), [], options);

[L, val, exitflag] = fmincon(@(L) objective_L(L, E, Pump, Signal), ...
    E.L, [], [], [], [], 0, 20, [], options);

function [y, dy] = objective(PdBm, E, Pump, Signal)
    % y = sum(Gain)
    % Gain gradient
    
    P = dBm2Watt(PdBm);
    Signal.P = P;
    [GaindB, ~, ~, dGaindB, dGaindB_L] = E.semi_analytical_gain(Pump, Signal);
    Gain = 10.^(GaindB/10);
          
%     y = -sum(10*log10(Gain));
    y = -sum(problem.step_approx(GaindB - problem.spanAttdB));
    
    dy = -sum(dGaindB.*problem.diff_step_approx(GaindB - problem.spanAttdB), 2);
    dy = log(10)/10*P.'.*dy;
end


function [y, dy] = objective_L(L, E, Pump, Signal)
    % y = sum(Gain)
    % Gain gradient
    
    E.L = L;
    [GaindB, ~, ~, dGaindB, dGaindB_L] = E.semi_analytical_gain(Pump, Signal);
     
    y = sum(GaindB);
    dy = sum(dGaindB_L);
%     dy = -sum(dGaindB.*problem.diff_step_approx(GaindB - problem.spanAttdB), 2);
end
end