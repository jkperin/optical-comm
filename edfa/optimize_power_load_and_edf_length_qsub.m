function C = optimize_power_load_and_edf_length_qsub(edf_type, totalPump, totalLengthKm, spanLengthKm)

addpath f/
addpath ../f/

verbose = false;

filename = sprintf('results/capacity_%s_totalPump=%sW_Ltotal=%s_Lspan=%skm.mat',...
        edf_type, totalPump, totalLengthKm, spanLengthKm);

filename = check_filename(filename) % verify if already exists and rename it if it does
    
% convert inputs to double (on cluster inputs are passed as strings)
if ~all(isnumeric([totalPump, totalLengthKm, spanLengthKm]))
    totalPump = str2double(totalPump);
    Ltotal = 1e3*str2double(totalLengthKm);
    Lspan = 1e3*str2double(spanLengthKm);
end

% EDF fiber
E = EDF(10, edf_type);

% Signal
df = 50e9;
dlamb = df2dlamb(df);
lamb = 1530e-9:dlamb:1565e-9;

Pon = 1e-4;
Signal = Channels(lamb, Pon, 'forward');
Pump = Channels(1480e-9, 30e-3, 'forward');

% SMF fiber
SMF = fiber(Lspan, @(lamb) 0.18, @(lamb) 0);
Namp = round(Ltotal/SMF.L); % number of amplifiers in chain
[~, spanAttdB] = SMF.link_attenuation(Signal.wavelength);

% Problem variables
problem.spanAttdB = spanAttdB;
problem.Namp = Namp;
problem.Pon = Pon;
problem.df = df;
        
% Pump power for each span
Pump.P = totalPump/Namp;

fprintf('Lspan = %d km, Pump.P = %.2f mW\n', Lspan/1e3, Pump.P*1e3);
        
% Optimize power load and EDF length
[E, Signal] = optimize_power_load_and_edf_length('particle swarm', E, Pump, Signal, problem, verbose);
                        
% Capacity calculation using numerical model
SignalOn = Channels(Signal.wavelength(Signal.P ~= 0), Signal.P(Signal.P ~= 0), 'forward');
ASEf = Channels(SignalOn.wavelength, 0, 'forward');
ASEb = Channels(SignalOn.wavelength, 0, 'backward');
GaindB_semi_analytical = E.semi_analytical_gain(Pump, SignalOn);
[GaindB, ~, ~, Pase, sol] = E.two_level_system(Pump, SignalOn, ASEf, ASEb, df, 50);
  
% Pase is ASE power in a bandwidth df
Gain = 10.^(GaindB/10); 
C = sum(log2(1 + SignalOn.P.*Gain./(Pase*(Namp-1))));

fprintf('Total spectrum efficiency = %.2f bits/s/Hz\n', C);

if verbose
    figure(205), hold on, box on
    plot(SignalOn.wavelength*1e9, GaindB, 'DisplayName', 'Numerical')
    plot(SignalOn.wavelength*1e9, GaindB_semi_analytical, 'DisplayName', 'Semi-analytical')
    xlabel('Wavelength (nm)')
    ylabel('Gain (dB)')
    legend('-DynamicLegend', 'Location', 'SouthEast')
end

% Save to file
save(filename)