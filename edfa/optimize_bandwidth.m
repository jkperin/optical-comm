%% Optimize capacity
clear, clc, close all

addpath ../f/

E = EDF(10, 'principles_type2');

df = 50e9;
dlamb = df2dlamb(df);
lamb = 1530e-9:dlamb:1565e-9;
L = 14e6; % total link distance in m
SMF = fiber(50e3, @(lamb) 0.18, @(lamb) 0);

Pon = 1e-4;
Signal = Channels(lamb, Pon, 'forward');
Pump = Channels(1480e-9, 30e-3, 'forward');

Ppump = (20:5:60)*1e-3;
Lspan = 30:10:80;
for l = 1:length(Lspan)
    SMF.L = Lspan(l)*1e3;
    N = L/SMF.L; % number of amplifiers in chain
    [~, SpanAttdB] = SMF.link_attenuation(Signal.wavelength);
    for p = 1:length(Ppump)
        fprintf('Lspan = %d km, Pump.P = %d mW\n', Lspan(l), Ppump(p)*1e3);
        Pump.P = Ppump(p);
        [Lopt(l, p), Signal] = optimize_edf_length(E, Pump, Signal, Pon, SpanAttdB, true);
        fprintf('- Optimal EDF length = %.2f\n', Lopt(l, p))
        
        % Capacity calculation using numerical model
        E.L = Lopt(l, p);
        SignalOn = Channels(Signal.wavelength(Signal.P ~= 0), Pon, 'forward');
        ASEf = Channels(SignalOn.wavelength, 0, 'forward');
        ASEb = Channels(SignalOn.wavelength, 0, 'backward');
        GaindB_semi_analytical = E.semi_analytical_gain(Pump, SignalOn);
        [GaindB, ~, ~, Pase, sol] = E.two_level_system(Pump, SignalOn, ASEf, ASEb, df, ceil(Lopt(l, p)/0.05));
        
        figure(204), hold on, box on
        plot(SignalOn.wavelength*1e9, GaindB, 'DisplayName', 'Numerical')
        plot(SignalOn.wavelength*1e9, GaindB_semi_analytical, 'DisplayName', 'Semi-analytical')
        xlabel('Wavelength (nm)')
        ylabel('Gain (dB)')
        legend('-DynamicLegend', 'Location', 'SouthEast')
        drawnow
        
        % Pase is ASE power in a bandwidth df
        Gain = 10.^(GaindB/10); 
        C(l, p) = sum(log2(1 + SignalOn.P.*Gain./(Pase*(N-1))));
    end
    
    figure(10), hold on, box on
    plot(Ppump*1e3, Lopt(l, :), '-o', 'DisplayName', sprintf('l = %d km', Lspan(l)))
    xlabel('Pump power (mW)')
    ylabel('Optimal EDF length (m)')
    legend('-dynamiclegend')
    
    figure(11), hold on, box on
    plot(Ppump*1e3, C(l, :), '-o', 'DisplayName', sprintf('l = %d km', Lspan(l)))
    xlabel('Pump power (mW)')
    ylabel('Maximum spectral efficiency (bit/s/Hz)')
    legend('-dynamiclegend')
    drawnow
end