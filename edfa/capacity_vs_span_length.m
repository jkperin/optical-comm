%% Optimize capacity
clear, clc, close all

addpath f/
addpath ../f/

E = EDF(10, 'principles_type3');

df = 50e9;
dlamb = df2dlamb(df);
lamb = 1530e-9:dlamb:1565e-9;
L = 14e6; % total link distance in m
SMF = fiber(50e3, @(lamb) 0.18, @(lamb) 0);

Pon = 1e-4;
Signal = Channels(lamb, Pon, 'forward');
Pump = Channels(1480e-9, 30e-3, 'forward');

problem.Pon = Pon;
problem.df = df;


TotalPump = 6:10;
% Ppump = (20:5:60)*1e-3;
Lspan = 30:5:80;
for p = 1:length(TotalPump)
    for l = 1:length(Lspan)
        SMF.L = Lspan(l)*1e3;
        N = L/SMF.L; % number of amplifiers in chain
        [~, spanAttdB] = SMF.link_attenuation(Signal.wavelength);
        problem.spanAttdB = spanAttdB;
        problem.Namp = N;
        
        % Pump power for each span
        Pump.P = TotalPump(p)/N;
        fprintf('Lspan = %d km, Pump.P = %.2f mW\n', Lspan(l), Pump.P*1e3);
        
        [E, Signal] = optimize_power_load_and_edf_length('particle swarm', E, Pump, Signal, problem, true);
        Lopt(p, l) = E.L;
                
        % Capacity calculation using numerical model
        E.L = Lopt(p, l);
        SignalOn = Channels(Signal.wavelength(Signal.P ~= 0), Signal.P(Signal.P ~= 0), 'forward');
        ASEf = Channels(SignalOn.wavelength, 0, 'forward');
        ASEb = Channels(SignalOn.wavelength, 0, 'backward');
        GaindB_semi_analytical = E.semi_analytical_gain(Pump, SignalOn);
        [GaindB, ~, ~, Pase, sol] = E.two_level_system(Pump, SignalOn, ASEf, ASEb, df, ceil(Lopt(p, l)/0.05));
        
        figure(204), hold on, box on
        plot(SignalOn.wavelength*1e9, GaindB, 'DisplayName', 'Numerical')
        plot(SignalOn.wavelength*1e9, GaindB_semi_analytical, 'DisplayName', 'Semi-analytical')
        xlabel('Wavelength (nm)')
        ylabel('Gain (dB)')
        legend('-DynamicLegend', 'Location', 'SouthEast')
        drawnow
        
        % Pase is ASE power in a bandwidth df
        Gain = 10.^(GaindB/10); 
        C(p, l) = sum(log2(1 + SignalOn.P.*Gain./(Pase*(N-1))));
    end
    
    figure(10), hold on, box on
    plot(Lspan, Lopt(p, :), '-o', 'DisplayName', sprintf('Total pump = %d W', TotalPump(p)))
    xlabel('Span length (km)')
    ylabel('Optimal EDF length (m)')
    legend('-dynamiclegend', 'Location', 'NorthWest')
    
    figure(11), hold on, box on
    plot(Lspan, C(p, :), '-o', 'DisplayName', sprintf('Total pump = %d W', TotalPump(p)))
    xlabel('Span length (km)')
    ylabel('Total spectral efficiency (bit/s/Hz)')
    legend('-dynamiclegend')
    drawnow
end