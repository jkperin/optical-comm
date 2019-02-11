%% Validate models with experimental data
clear, clc, close all

addpath data/
addpath f/
addpath ../f/

experiment = 'Pin=-13dBm'; % which power values to load
experiment_April_6_18 % load appropriate power values

E = EDF(8.64, 'corning_exp');
E.core_radius = 1.19e-6;
E.doping_radius = 0.98e-6; % Er3+ core radius. e.g., 1.2 um in [1, Table 1], 1.05um in [4, pg 156]
E.rho0 = 6.68e18; % Er3+ concentraction (cm^3), e.g., 0.7e19 in [4, pg 156]
E.NA = 0.23; % numerical aperture, e.g., 0.28 in [4, pg. 156]
E.tau = 10e-3; % metastable lifetime in s        
E.alphap_980nm = 4.272;
E.gp_980nm = 0;
nsp_correction = 1.4;

E.excess_loss = 0.28; % -13 dBm
% E.excess_loss = 0.31; % -16 dBm
% E.excess_loss = 0.35; % -19 dBm
offset = 0;

% E.excess_loss = 0.35;
% offset = 2.8;
% coupling_loss = 0.2;


Signal = Channels(wavelength.'*1e-9, dBm2Watt(Pin.' + offset), 'forward');

leg = {};
for k = 1:length(PpumpmW)
    Pump = Channels(980e-9, PpumpmW(k)*1e-3, 'forward');
    ASEf = Channels(Signal.wavelength, 0, 'forward');
    ASEb = Channels(Signal.wavelength, 0, 'backward');  

%     [x,fval,exitflag] = fminunc(@(x) fit_excess_loss(x, E, Pump, Signal, experimental_GaindB(:, k), experimental_PasedBm(:, k) + offset), [0.3 0.3 0 0])
%     [x,fval,exitflag] = fminunc(@(x) fit_EDF_length(x, E, Pump, wavelength, Pin, BWref,...
%         experimental_GaindB(:, k), experimental_PasedBm(:, k)), [9 0.2 offset 0])
%     
%     E.excess_loss = x(1);
%     pump_coupling_loss_dB = x(2);
%     signal_coupling_loss_dB = x(3);
%     ASE_coupling_loss_dB = x(4);
%     
%     Pump.P = dBm2Watt(Watt2dBm(PpumpmW(k)*1e-3) - pump_coupling_loss_dB);
%     Signal.P = dBm2Watt(Pin.' + offset - signal_coupling_loss_dB);
%     experimental_PasedBm(:, k) = experimental_PasedBm(:, k) - ASE_coupling_loss_dB;
     
%     [GaindB, Pump_out, Psignal_out, Pase, sol] = E.SHB(Pump, Signal, ASEf, ASEb, BWref, 'three-level', 100, false);
    [GaindB, Pump_out, Psignal_out, Pase, sol] = E.propagate(Pump, Signal, ASEf, ASEb, BWref, 'three-level', 50, false);
    GaindB_sa = E.semi_analytical_gain(Pump, Signal);
    Pase_sa = E.analytical_ASE_PSD(Pump, Signal, nsp_correction);
    Pase_sa = Pase_sa*BWref;
    
    [~, NFdB] = E.excess_noise(Pump, Signal, false);

    figure(1), hold on, box on
    hplot1(k) = plot(wavelength, GaindB, '-');
    plot(wavelength, GaindB_sa, '--', 'Color', get(hplot1(k), 'Color'))
    plot(wavelength, experimental_GaindB(:, k), '.', 'Color', get(hplot1(k), 'Color'))
    xlabel('Wavelength (nm)', 'FontSize', 12)
    ylabel('Gain (dB)', 'FontSize', 12)
    set(gca, 'FontSize', 12)
    
%     figure(11), hold on, box on
%     plot(wavelength, GaindB - experimental_GaindB(:, k).', '-');
%     xlabel('Wavelength (nm)', 'FontSize', 12)
%     ylabel('Gain difference (dB)', 'FontSize', 12)
%     set(gca, 'FontSize', 12)
     
%     figure(10), hold on, box on
%     hplot1(k) = plot(wavelength, GaindB, '-');
    experimental_GaindB_unbiased = 10*log10(10.^(experimental_GaindB(:, k)/10) - dBm2Watt(experimental_PasedBm(:, k) + offset)./dBm2Watt(Pin+offset));
%     plot(wavelength, experimental_GaindB_unbiased, '.', 'Color', get(hplot1(k), 'Color'))
%     xlabel('Wavelength (nm)', 'FontSize', 12)
%     ylabel('Gain (dB)', 'FontSize', 12)
%     set(gca, 'FontSize', 12)

    figure(2), hold on, box on
    hplot2(k) = plot(wavelength, Watt2dBm(Pase), '-');
%     plot(wavelength, Watt2dBm(Pase_sa), '--', 'Color', get(hplot2(k), 'Color'))
    plot(wavelength, experimental_PasedBm(:, k) + offset, '.', 'Color', get(hplot2(k), 'Color'))
    xlabel('Wavelength (nm)', 'FontSize', 12)
    ylabel('ASE power (dBm)', 'FontSize', 12)
    set(gca, 'FontSize', 12)
%     
%     figure(3), hold on, box on
%     hplot3(k) = plot(wavelength, NFdB, '-');
%     NFdB_exp = 10*log10(dBm2Watt(experimental_PasedBm(:, k).' + offset)./(BWref*Signal.Ephoton.*10.^(experimental_GaindB_unbiased.'/10)));
%     plot(wavelength, NFdB_exp, '.', 'Color', get(hplot3(k), 'Color'))
%     xlabel('Wavelength (nm)', 'FontSize', 12)
%     ylabel('Noise figure (dB)', 'FontSize', 12)
%     set(gca, 'FontSize', 12)

%     figure(3), hold on, box on
%     hplot3(k) = plot(wavelength, Watt2dBm(Psignal_out) - Watt2dBm(Pase), '-');
%     plot(wavelength, Pin + experimental_GaindB(:, k) - experimental_PasedBm(:, k), '.', 'Color', get(hplot3(k), 'Color'))
%     xlabel('Wavelength (nm)')
%     ylabel('OSNR 0.1nm (dB)')
%     drawnow
    leg = [leg sprintf('Pump = %.2f mW', PpumpmW(k))];

    drawnow
end

figure(1), legend(hplot1, leg, 'Location', 'SouthEast');
% ylim([2 20])
% m = matlab2tikz(gca);
% m.write_tables('gain_experimental', 'same x')
figure(2), legend(hplot2, leg, 'Location', 'SouthEast');
% ylim([-50 -34])
% figure(3), legend(hplot3, leg, 'Location', 'NorthEast');
% ylim([34 39])
% m = matlab2tikz(gca);
% m.write_tables('ase_experimental', 'same x')

