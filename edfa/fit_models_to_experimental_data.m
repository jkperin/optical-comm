%% Fit models to experimental data
clear, clc, close all

addpath data/
addpath f/
addpath ../f/

%% Import experimental data from file
data_set_name = '-10.5dBm input per ch';

Experiment = load_experimental_data(data_set_name);

data = Experiment(data_set_name); % select desired experiment

%% 

E = EDF(6, 'corning high na');

% Parameters that minimize mean square error of Gain and ASE for each data set
% Computed using function f/fit_param_to_data.m
% -16.5 dBm
% E.excess_loss = 0.4114; 
% nsp_correction = 1.30006;

% -14.5 dBm
% E.excess_loss = 0.3907; 
% nsp_correction = 1.3987;

% -12.5 dBm
% E.excess_loss = 0.3768; 
% nsp_correction = 1.5488;

% -10.5 dBm
% E.excess_loss = 0.4028; 
% nsp_correction = 1.81557;

% Selected
E.excess_loss = 0.4011; 
nsp_correction = 1.3494;


ReferenceBandwidth = 12.5e9; % Bandwidth for noise power measurement

wavelength = data(1).wavelengthnm.';
input_power_mW = dBm2Watt(data(1).InputPowerdBm.');

Signal = Channels(wavelength * 1e-9, input_power_mW, 'forward');

leg = {};
for k = 1:length(data)
    pump_power_mW = data(k).PumpPowermW;
    
    Pump = Channels(980e-9, pump_power_mW * 1e-3, 'forward');
    ASEf = Channels(Signal.wavelength, 0, 'forward');
    ASEb = Channels(Signal.wavelength, 0, 'backward');  

    [GaindB, Pump_out, Psignal_out, Pase, sol] = E.propagate(Pump, Signal,...
        ASEf, ASEb, ReferenceBandwidth, 'two-level', 50, false);
    GaindB_sa = E.semi_analytical_gain(Pump, Signal);
    Pase_sa = E.analytical_ASE_PSD(Pump, Signal, nsp_correction);
    Pase_sa = Pase_sa*ReferenceBandwidth;

    experimental_Gain = 10.^(data(k).GaindB.'/10);
    unbiased_experimental_GaindB = 10*log10(experimental_Gain - dBm2Watt(data(k).ASEdBm.')./Signal.P);
    
    figure(1), hold on, box on
    hplot1(k) = plot(wavelength, GaindB, '-');
%     plot(wavelength, GaindB_sa, '--', 'Color', get(hplot1(k), 'Color'))
    plot(wavelength, unbiased_experimental_GaindB, '.', 'Color', get(hplot1(k), 'Color'))
    xlabel('Wavelength (nm)', 'FontSize', 12)
    ylabel('Gain (dB)', 'FontSize', 12)
    set(gca, 'FontSize', 12)
    
    figure(2), hold on, box on
    hplot2(k) = plot(wavelength, Watt2dBm(Pase), '-');
%     plot(wavelength, Watt2dBm(Pase_sa), '--', 'Color', get(hplot2(k), 'Color'))
    plot(wavelength, data(k).ASEdBm, '.', 'Color', get(hplot2(k), 'Color'))
    xlabel('Wavelength (nm)', 'FontSize', 12)
    ylabel('ASE power (dBm)', 'FontSize', 12)
    set(gca, 'FontSize', 12)

    leg = [leg sprintf('Pump = %d mW', data(k).NominalPumpPowermW)];

    drawnow
end

figure(1), legend(hplot1, leg, 'Location', 'NorthEast');
% ylim([5 20])
% m = matlab2tikz(gca);
% m.write_tables('gain_experimental', 'same x')
figure(2), legend(hplot2, leg, 'Location', 'NorthEast');
% ylim([-50 -34])
% figure(3), legend(hplot3, leg, 'Location', 'NorthEast');
% ylim([34 39])
% m = matlab2tikz(gca);
% m.write_tables('ase_experimental', 'same x')

