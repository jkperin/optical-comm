%% Model comparisons
clear, clc, close all

addpath data/
addpath f/
addpath ../f/

E = EDF(8, 'corning_type1');
% E = EDF(10, 'giles_ge:silicate');

PonVec = -20:-10;
excess_noise_correction = 1.4;

for k = 1:length(PonVec)
    Pon = dBm2Watt(PonVec(k));
    Pump = Channels(980e-9, 60e-3, 'forward'); %
    Signal = Channels(linspace(1530, 1565, 88)*1e-9, Pon, 'forward');
    ASEf = Channels(Signal.wavelength, 0, 'forward');
    ASEb = Channels(Signal.wavelength, 0, 'backward');

    dlamb = Signal.wavelength(2)-Signal.wavelength(1);
    df = E.c/Signal.wavelength(1)-E.c/Signal.wavelength(2);
%     df = 50e9;
        
    GaindB_semi_analytical = E.semi_analytical_gain(Pump, Signal);
    Pase_analytical = E.analytical_ASE_PSD(Pump, Signal, excess_noise_correction)*df; % ASE power

    [GaindB, Ppump_out, Psignal_out, Pase, sol] = E.two_level_system(Pump, Signal, ASEf, ASEb, df, 50, false);
        
    figure(1), hold on, box on
    hplot = plot(Signal.wavelength*1e9, GaindB);
    plot(Signal.wavelength*1e9, GaindB_semi_analytical, '--', 'Color', get(hplot, 'Color'))
    xlabel('Wavelength (nm)', 'FontSize', 12)
    ylabel('Gain (dB)', 'FontSize', 12)
        
    figure(2), hold on, box on
    hplot = plot(Signal.wavelength*1e9, Watt2dBm(Pase));
    plot(Signal.wavelength*1e9, Watt2dBm(Pase_analytical), '--', 'Color', get(hplot, 'Color'))
    xlabel('Wavelength (nm)', 'FontSize', 12)
    ylabel('ASE (dBm)', 'FontSize', 12)
    
    figure(3), hold on, box on
    hplot = plot(Signal.wavelength*1e9, Watt2dBm(Psignal_out));
    plot(Signal.wavelength*1e9, Signal.PdBm + GaindB_semi_analytical, '--', 'Color', get(hplot, 'Color'))
    xlabel('Wavelength (nm)', 'FontSize', 12)
    ylabel('Output power (dBm)', 'FontSize', 12)
    drawnow
end
