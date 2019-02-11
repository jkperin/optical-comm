%% Animation of particle swarm optimization

addpath data/
addpath f/
addpath ../f/

% EDF fiber
E = EDF(10, 'corning_type1');
E.excess_loss = 0;

% Pump & Signal
nonlinear_coeff_file = sprintf('../f/GN_model_coeff_spanLengthkm=%dkm_Df=%dGHz.mat', 50, 50);
NCOEFF = load(nonlinear_coeff_file);
lamb = NCOEFF.lamb;
    
Signal = Channels(lamb, 0, 'forward');
Pump = Channels(980e-9, 60e-3, 'forward');

% SMF fiber
SMF = fiber(50e3, @(l) 0.165*ones(size(l)), @(l) 20.4e-6*ones(size(l)));
SMF.gamma = 0.8e-3;
[~, spanAttdB] = SMF.link_attenuation(1550e-9); % same attenuation for all wavelengths
spanAttdB = spanAttdB + 1.5; % adds 1.5 dB of margin

% Problem variables
problem.spanAttdB = spanAttdB;
problem.df = 50e9;
problem.Gap = 10^(-1/10);
problem.Namp = 287;
problem.step_approx = @(x) 0.5*(tanh(2*x) + 1); % Smoothing factor = 2
problem.diff_step_approx = @(x) sech(2*x).^2; % first derivative (used for computing gradient)
problem.excess_noise_correction = 1.55; 
problem.SwarmSize = min(100, 20*(Signal.N+1));
problem.nonlinearity = false;
problem.nonlinear_coeff = NCOEFF.nonlinear_coeff;
problem.epsilon = 0.07; % From Fig. 17 of P. Poggiolini and I. Paper, “The GN Model
% of Non-Linear Propagation in Uncompensated Coherent Optical Systems,” 
% J. Lightw. Technol., vol. 30, no. 24, pp. 3857–3879, 2012.

% Leff = SMF.effective_length(1550e-9);
% BW = 100*problem.df; % about 100 channels
% problem.epsilon = 0.3*log(1 + 6*Leff/(SMF.L*asinh(pi^2/2*abs(SMF.beta2(1550e-9))*Leff*BW^2))); % From Eq. 23 of P. Poggiolini and I. Paper, “The GN Model
% % of Non-Linear Propagation in Uncompensated Coherent Optical Systems,” 
% % J. Lightw. Technol., vol. 30, no. 24, pp. 3857–3879, 2012.


% Define power cap according to conservation of energy and add 3 dB as
% safety margin (multiplication by 2)
problem.Pon = 2*(1/Signal.N)*Pump.wavelength*Pump.P/(min(Signal.wavelength)*(10^(spanAttdB/10)-1));
Signal.P(1:end) = problem.Pon;

fprintf('Power cap = %.5f mW (%.3f dBm)\n', 1e3*problem.Pon, Watt2dBm(problem.Pon));

problem.excess_noise = E.analytical_excess_noise(Pump, Signal);
problem.excess_noise = problem.excess_noise*problem.excess_noise_correction;

%% Optimize EDF length and power allocation jointly
SwarmSize = 20;
SwarmMatrix = zeros(SwarmSize, Signal.N);
for k = 1:SwarmSize
    SwarmMatrix(k, :) = (Watt2dBm(problem.Pon) - k/2)*ones(1, Signal.N);
end

la = -Inf*ones(1, Signal.N); % lower bound
lb = Watt2dBm(problem.Pon)*ones(1, Signal.N); % upper bound
options = optimoptions('particleswarm', 'Display', 'iter', 'UseParallel', true,...
        'MaxStallTime', 60, 'MaxStallIterations', 100, 'SwarmSize', SwarmSize,...
        'InitialSwarmSpan', 10*ones(1, Signal.N),...
        'InitialSwarmMatrix', SwarmMatrix,...
        'PlotFcn', @(x,y) pso_animation_plot_function(x, y, Signal.lnm)); % SwarmSpan for fiber length is 20m and 10 dB for signal power

[X, relaxed_SE, exitflag] = particleswarm(@(X) -capacity_linear_regime_relaxed([6.6 X], E, Pump, Signal, problem),...
    Signal.N, la, lb, options);

