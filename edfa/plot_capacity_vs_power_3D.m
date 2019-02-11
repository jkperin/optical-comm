%% Calculate capacity for small variation in power of two channels and plot the results 
clear, clc, close all

addpath ../
addpath ../f/
addpath ../../f/
addpath ../../f/DERIVEST

S = load('../results/capacity_vs_pump_power_EDF=corning_type1_pump=980nm_L=286_x_50km/capacity_vs_pump_power_EDF=corning_type1_pump=150mW_980nm_L=286_x_50km.mat');

E = S.nlin.E;
Signal = S.nlin.S;
Pump = S.Pump;
problem = S.problem;
problem.step_approx = @(x) 0.5*(tanh(2*x) + 1); % Smoothing factor = 2
problem.diff_step_approx = @(x) sech(2*x).^2; % First derivative of step_approx
problem.excess_noise = E.analytical_excess_noise(Pump, Signal);
problem.excess_noise = problem.excess_noise*problem.excess_noise_correction;

N = Signal.N;
Signal = Signal.sample(Signal.P ~= 0);
Ndiscard = N - Signal.N;
problem.nonlinear_coeff{1} = problem.nonlinear_coeff{1}(Ndiscard+1:end-Ndiscard, Ndiscard+1:end-Ndiscard);
problem.nonlinear_coeff{2} = problem.nonlinear_coeff{2}(Ndiscard+1:end-Ndiscard, Ndiscard+1:end-Ndiscard);
problem.nonlinear_coeff{3} = problem.nonlinear_coeff{3}(Ndiscard+1:end-Ndiscard, Ndiscard+1:end-Ndiscard);
problem.excess_noise = E.analytical_excess_noise(Pump, Signal);
problem.excess_noise = problem.excess_noise*problem.excess_noise_correction;

% 
ch_idx = [40 81]; % index of channels to be varied
Npoints = 50;
deltaPdBm = 10;

X = linspace(-deltaPdBm/2, deltaPdBm/2, Npoints);
Y = linspace(-deltaPdBm/2, deltaPdBm/2, Npoints);
[X, Y] = meshgrid(X, Y);

PdBm = Signal.PdBm;

SE = zeros(size(X));
SE2 = zeros(size(X));
for i = 1:Npoints
    fprintf('(%d, %d)\n', i, j)
    for j = 1:Npoints
        PdBm(ch_idx) = Signal.PdBm(ch_idx) + [X(i, j) Y(i, j)];
        
        SE(i, j) = capacity_nonlinear_regime_relaxed([E.L PdBm], E, Pump, Signal, problem);
        SE2(i, j) = capacity_nonlinear_regime_relaxed_all_on(PdBm, Signal, problem);
    end
end

figure, 
mesh(X, Y, SE)

figure, 
contour(X, Y, SE)
    
figure, 
mesh(X, Y, SE2)

figure, 
contour(X, Y, SE2)



