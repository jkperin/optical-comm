function capacity_vs_pump_power_slurm(task)
%% Optimize channel power and EDF length to maximize capacity for a given pump power
% task is an integer (passed as string) that indexes the array taskList

addpath data/
addpath f/
addpath ../f/

verbose = false;

% Select task
taskList = [30:5:150 200:100:500]; % Variable to be modified in different system calls
pumpPowermW = taskList(round(str2double(task)));
pumpPower = 1e-3*pumpPowermW;

% Other input parameters
Nspans = 286;
spanLengthKm = 50;

% EDF fiber
E = EDF(10, 'corning_type1');

% Pump & Signal
df = 50e9;
dlamb = df2dlamb(df);
lamb = 1522e-9:dlamb:1582e-9;
Signal = Channels(lamb, 0, 'forward');
Pump = Channels(980e-9, pumpPower, 'forward');

% SMF fiber
SMF = fiber(spanLengthKm*1e3, @(lamb) 0.18, @(lamb) 0);
SMF.gamma = 0.8e-3;
[~, spanAttdB] = SMF.link_attenuation(Signal.wavelength);

% Filename
filename = sprintf('results/capacity_vs_pump_power_EDF=%s_pump=%dmW_%dnm_L=%d_x_%dkm.mat',...
        E.type, pumpPowermW, Pump.wavelength*1e9, Nspans, spanLengthKm);
filename = check_filename(filename); % verify if already exists and rename it if it does
disp(filename) 

% Problem variables
problem.spanAttdB = spanAttdB;
problem.df = df;
problem.Namp = Nspans;
problem.step_approx = @(x) 0.5*(tanh(2*x) + 1); % Smoothing factor = 2
problem.diff_step_approx = @(x) sech(2*x).^2; % first derivative (used for computing gradient)
problem.excess_noise_correction = 1.4; % 1.2 for 980nm, 1.6 for 1480nm
problem.SwarmSize = min(100, 20*(Signal.N+1));
problem.nonlinearity = true;
S = load('../f/GN_model_coeff_spanLengthkm=50km_Df=50GHz.mat');
problem.nonlinear_coeff = S.nonlinear_coeff;
problem.epsilon = 0.05; % From Fig. 17 of P. Poggiolini and I. Paper, “The GN Model
% of Non-Linear Propagation in Uncompensated Coherent Optical Systems,” 
% J. Lightw. Technol., vol. 30, no. 24, pp. 3857–3879, 2012.

% Define power cap according to conservation of energy and add 3 dB as
% safety margin (multiplication by 2)
problem.Pon = 2*(1/Signal.N)*Pump.wavelength*Pump.P/(min(Signal.wavelength)*(10^(spanAttdB/10)-1));
Signal.P(1:end) = problem.Pon;

fprintf('Power cap = %.5f mW (%.3f dBm)\n', 1e3*problem.Pon, Watt2dBm(problem.Pon));

%% Optimize power load and EDF length 
try    
    %% Linear regime
    disp('== Linear regime')
    problem.nonlinearity = false;
    [lin.E, lin.S, lin.exitflag, lin.num, lin.approx] = ...
        optimize_power_load_and_edf_length('particle swarm', E, Pump, Signal, problem, verbose);
    
    failed_lin = lin;   
    if lin.exitflag ~= 1 % try one more time
        fprintf('Optimization ended with exitflag = %d\n. Trying one more time...\n', lin.exitflag)
        [lin.E, lin.S, lin.exitflag, lin.num, lin.approx] = ...
            optimize_power_load_and_edf_length('particle swarm', E, Pump, Signal, problem, verbose);
    end
catch e
    warning(e.message)
end
  
try
    %% Nonlinear regime
    problem.nonlinearity = true;
    disp('== Nonlinear regime (global optimization)')
    [nlin.E, nlin.S, nlin.exitflag, nlin.num, nlin.approx] = ...
        optimize_power_load_and_edf_length('particle swarm', E, Pump, Signal, problem, verbose);

    % local optimization after particle swarm is done
    problem.nonlinearity = true;
    disp('== Nonlinear regime (local optimization)')
    [nlin_local.E, nlin_local.S, nlin_local.exitflag, nlin_local.num, nlin_local.approx] = ...
        optimize_power_load_and_edf_length('local', nlin.E, Pump, nlin.S, problem, verbose);
    
catch e
    warning(e.message)
end

% Save to file
save(filename)