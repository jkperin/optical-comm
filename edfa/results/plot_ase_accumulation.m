clc, clear, close all

addpath ../../f/

% Paper simulations
% addpath capacity_vs_pump_power_highNA_final
% load('sim_ase_acc_Ppump=60mW_GFF_period=1_partial_GFF=0.mat')

addpath ase_accumulation
% load('sim_ase_acc_Ppump=60mW_GFF_period=1_partial_GFF=0')
load('sim_ase_acc_Ppump=60mW_GFF_period=3_partial_GFF=0_SpanLoss-2dB')

% idx = SignalIn.P ~= OffPower;
% idx = 1:SignalIn.N;
idx = [1 100 200 220];

figure(1), hold on, box on
plot(SignalIn.lnm, nGaindB(:, idx), 'LineWidth', 2)
xlabel('Wavelength (nm)', 'FontSize', 12)
ylabel('Gain (dB)', 'FontSize', 12)
legend('Span = 1', 'Span = 100', 'Span = 200', 'Span = 220')
set(gca, 'FontSize', 12)
% m = matlab2tikz(gca);
% m.write_tables('sim_gain', 'same x')

% figure(2), hold on, box on
% plot(SignalIn.lnm, nGFF(:, idx), 'LineWidth', 2)
% xlabel('Wavelength (nm)', 'FontSize', 12)
% ylabel('Ideal GFF (dB)', 'FontSize', 12)
% legend('Span = 1', 'Span = 100', 'Span = 200', 'Span = 220')
% set(gca, 'FontSize', 12)
% m = matlab2tikz(gca);
% m.write_tables('sim_gff', 'same x')

% 
figure(3), hold on, box on
plot(SignalIn.lnm, nASE(:, idx), 'LineWidth', 2)
xlabel('Wavelength (nm)', 'FontSize', 12)
ylabel('Accumulated ASE (dBm)', 'FontSize', 12)
legend('Span = 1', 'Span = 100', 'Span = 200', 'Span = 220')
set(gca, 'FontSize', 12)
% m = matlab2tikz(gca);
% m.write_tables('sim_ase', 'same x')

figure(4), hold on, box on
plot(SignalIn.lnm, SE, 'LineWidth', 2)
plot(SignalIn.lnm, S.lin.num.SE, 'LineWidth', 2)
xlabel('Wavelength (nm)', 'FontSize', 12)
ylabel('Spectral efficiency (bits/s/Hz)', 'FontSize', 12)
set(gca, 'FontSize', 12)
legend('Simulation', 'Optimization')
% m = matlab2tikz(gca);
% m.write_tables('sim_se', 'same x')