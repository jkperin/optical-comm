%% TE subcomm experiment
clear, clc

addpath data/
addpath f/
addpath ../f/

E = EDF(7, 'corning_type1');

df = 33.3e9;  
dlamb = df2dlamb(df);
lamb = 1525e-9:dlamb:1570e-9;
L = 14.35e6;
SMF = fiber(46e3, @(lamb) 0.211, @(lamb) 0);
SMF.gamma = 0.8e-3;
Namp = round(L/SMF.L);

Signal = Channels(lamb, 0, 'forward');
Pump = Channels(980e-9, 60e-3, 'forward');

[~, spanAttdB] = SMF.link_attenuation(Signal.wavelength);

problem.Pon = 2*(1/Signal.N)*Pump.wavelength*Pump.P/(min(Signal.wavelength)*(10^(spanAttdB/10)-1));
problem.spanAttdB = spanAttdB;
problem.df = df;
problem.Gap = 10^(-1/10);
problem.Namp = Namp;
problem.step_approx = @(x) 0.5*(tanh(2*x) + 1); % Smoothing factor = 2
problem.diff_step_approx = @(x) sech(2*x).^2; % first derivative (used for computing gradient)
problem.excess_noise_correction = 1.4; % 1.2 for 980nm, 1.6 for 1480nm
problem.SwarmSize = min(200, 20*(Signal.N+1));
problem.nonlinearity = true;
S = load('../f/GN_model_coeff_spanLengthkm=50km_Df=33GHz');
problem.nonlinear_coeff = S.nonlinear_coeff;
problem.epsilon = 0.05; % From Fig. 17 of P. Poggiolini and I. Paper, “The GN Model
% of Non-Linear Propagation in Uncompensated Coherent Optical Systems,” 
% J. Lightw. Technol., vol. 30, no. 24, pp. 3857–3879, 2012.

%% TE subcomm experiment
% Output power = 12 dBm for 82 channels
% Measured OSNR was 11 dB in 0.1 nm
% Span attenuation = 9.7, so gain is about 10 dB in average
Signal.P(Signal.wavelength > 1539e-9 & Signal.wavelength < 1561e-9) = dBm2Watt(12)/(82*10);

% Flat power allocation
[E, Signal, ~, num, approx] ... 
    = optimize_power_load_and_edf_length('none', E, Pump, Signal, problem, true);

% Power optimization
% [Eopt_pswarm, SignalOn_pswarm, exitflag_pswarm, num_pswarm, approx_pswarm] ... 
%     = optimize_power_load_and_edf_length('particle swarm', E, Pump, Signal, problem, false);
% 
% [Eopt_local, SignalOn_local, exitflag_local, num_local, approx_local] ... 
%     = optimize_power_load_and_edf_length('local', Eopt_pswarm, Pump, SignalOn_pswarm, problem, true);

