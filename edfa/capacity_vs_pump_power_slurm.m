function capacity_vs_pump_power_slurm(spanLengthKm, Nspans, task)
%% Optimize channel power and EDF length to maximize capacity for a given pump power
% task is an integer (passed as string) that indexes the array taskList

addpath data/
addpath f/
addpath ../f/

verbose = false;

% Select task
taskList = [20:5:200 225:25:400]; % Variable to be modified in different system calls
pumpPowermW = taskList(round(str2double(task)));
pumpPower = 1e-3*pumpPowermW;

% Other input parameters
spanLengthKm = round(str2double(spanLengthKm));
Nspans = round(str2double(Nspans));
freqSpacingGHz = 50;

% EDF fiber
E = EDF(10, 'corning high NA');

% Pump & Signal
nonlinear_coeff_file = sprintf('../f/GN_model_coeff_spanLengthkm=%dkm_Df=%dGHz.mat', spanLengthKm, freqSpacingGHz);
NCOEFF = load(nonlinear_coeff_file);
lamb = NCOEFF.lamb;
lamb = lamb(lamb <= 1570e-9); % limit wavelengths according to fiber data

Signal = Channels(lamb, 0, 'forward');
Pump = Channels(980e-9, pumpPower, 'forward');

% SMF fiber
SMF = fiber(spanLengthKm*1e3, @(l) 0.165*ones(size(l)), @(l) 20.4e-6*ones(size(l)));
SMF.gamma = 0.8e-3;
[~, spanAttdB] = SMF.link_attenuation(1550e-9); % same attenuation for all wavelengths
spanAttdB = spanAttdB + 1.5; % adds 1.5 dB of margin

% Filename
filename = sprintf('results/capacity_vs_pump_power_EDF=%s_pump=%dmW_%dnm_ChDf=%dGHz_L=%d_x_%dkm.mat',...
        E.type, pumpPowermW, round(Pump.wavelength*1e9), freqSpacingGHz, Nspans, spanLengthKm);
filename = check_filename(filename); % verify if already exists and rename it if it does
disp(filename) 

% Problem variables
problem.spanAttdB = spanAttdB;
problem.df = freqSpacingGHz*1e9;
problem.Gap = 10^(-1/10);
problem.Namp = Nspans;
problem.step_approx = @(x) 0.5*(tanh(2*x) + 1); % Smoothing factor = 2
problem.diff_step_approx = @(x) sech(2*x).^2; % first derivative (used for computing gradient)
problem.excess_noise_correction = 1.55; 
problem.SwarmSize = min(300, 20*(Signal.N+1));
problem.nonlinearity = true;
problem.nonlinear_coeff = NCOEFF.nonlinear_coeff;
problem.epsilon = 0.07; % From Fig. 17 of P. Poggiolini and I. Paper, “The GN Model
% of Non-Linear Propagation in Uncompensated Coherent Optical Systems,” 
% J. Lightw. Technol., vol. 30, no. 24, pp. 3857–3879, 2012.

% Leff = SMF.effective_length(1550e-9);
% BW = 100*problem.df; % about 100 channels
% problem.epsilon = 0.3*log(1 + 6*Leff/(SMF.L*asinh(pi^2/2*abs(SMF.beta2(1550e-9))*Leff*BW^2))); % From Eq. 23 of P. Poggiolini and I. Paper, “The GN Model
% % of Non-Linear Propagation in Uncompensated Coherent Optical Systems,” 
% % J. Lightw. Technol., vol. 30, no. 24, pp. 3857–3879, 2012.

if freqSpacingGHz == 50
    options.AdaptationConstant = 0.1; 
else
    options.AdaptationConstant = 0.01; 
end
options.FiniteDiffStepSize = 1e-6;
options.MaxIterations = 50;
options.AbsTol = 1e-3;
options.MinStep = 1e-4;
problem.saddle_free_newton.options = options;

% Define power cap according to conservation of energy and add 3 dB as
% safety margin (multiplication by 2)
problem.Pon = (1/(0.8*Signal.N))*Pump.wavelength*Pump.P/(min(Signal.wavelength)*(10^(spanAttdB/10)-1));
Signal.P(1:end) = problem.Pon;

fprintf('Power cap = %.5f mW (%.3f dBm)\n', 1e3*problem.Pon, Watt2dBm(problem.Pon));

%% Optimize power load and EDF length 
try    
    %% Linear regime
    disp('== Linear regime')
    problem.nonlinearity = false;
    [lin.E, lin.S, lin.exitflag, lin.num, lin.approx] = ...
        optimize_power_load_and_edf_length('particle swarm', E, Pump, Signal, problem, verbose);
    
    if lin.exitflag ~= 1 % try one more time
        failed_lin = lin;   
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

    try
        % local optimization after particle swarm is done
        disp('== Nonlinear regime (unconstrained optimization)')
        [nlin_unc.E, nlin_unc.S, nlin_unc.exitflag, nlin_unc.num, nlin_unc.approx] = ...
            optimize_power_load_and_edf_length('local unconstrained', nlin.E, Pump, nlin.S, problem, verbose);
    catch ee
        warning(ee.message)
    end

    try
        % Continue optimization using Saddle-free Newton method
        disp('== Nonlinear regime (Saddle-free Newton)')
        [nlin_sfn.E, nlin_sfn.S, nlin_sfn.exitflag, nlin_sfn.num, nlin_sfn.approx] = ...
            optimize_power_load_and_edf_length('saddle-free newton', nlin.E, Pump, nlin.S, problem, verbose);
    catch ee
        warning(ee.message)
    end
    
catch e
    warning(e.message)
end

% Save to file
save(filename)
