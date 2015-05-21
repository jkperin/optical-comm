clear, clc, close all

%% Evaluate different noise regimes
% PSDs are unilateral, unless otherwise indicated

PrxdBm = linspace(-16,8);
Prx = 1e-3*10.^(PrxdBm/10);


% Simulation parameters
Nc = 64;
Nu = 52;
rclip = 3.7;
fs = [65.8e9 43.9e9];  % sampling rate
ENOB = 6;

ros = Nc/Nu;

BW = fs/ros; % two-sided noise bandwidth BW = fs/ros

%% Thermal noise
NEP = 30e-12; % A/sqrt(Hz)
N0 = NEP^2;

Sth = @(P) 10*log10(NEP^2)*ones(size(P));

%% Shot Noise
R = 1;
q = 1.60217657e-19;      % electron charge (C)
Id = 0;                  % dark current

Sshot = @(P) 2*q*(R*P + Id);     % one-sided shot noise PSD

%% Intensity noise
RIN = -150;     % dB/Hz
% RIN must be in dB/Hz
Srin = @(RIN, P) 10^(RIN/10)*2*R^2*P.^2;

%
% Pn = (Prx/rclip).^2/Nu;

varQ = Prx.^2/(3*2^(2*ENOB));

SNR16QAM = 10*log10(ros*(Prx/rclip).^2./(fs(1)/2*(N0 + 2*q*R*Prx + Srin(RIN, Prx))));
SNR64QAM = 10*log10(ros*(Prx/rclip).^2./(fs(2)/2*(N0 + 2*q*R*Prx + Srin(RIN, Prx))));

SNR16QAMQ = 10*log10(ros*(Prx/rclip).^2./(fs(1)/2*(N0 + 2*q*R*Prx + Srin(RIN, Prx)) + 2*varQ));
SNR64QAMQ = 10*log10(ros*(Prx/rclip).^2./(fs(2)/2*(N0 + 2*q*R*Prx + Srin(RIN, Prx)) + 2*varQ));

% SNR16QAM = 10*log10(Pn./(1/Nc*BW(1)*(N0/2 + Sshot(Prx)/2 + Srin(RIN, Prx)/2)));
% SNR64QAM = 10*log10(Pn./(1/Nc*BW(2)*(N0/2 + Sshot(Prx)/2 + Srin(RIN, Prx)/2)));
% 
% SNR16QAMQ = 10*log10(Pn./(1/Nc*BW(1)*(N0/2 + Sshot(Prx)/2 + Srin(RIN, Prx)/2) + 1/Nc*2*varQ));
% SNR64QAMQ = 10*log10(Pn./(1/Nc*BW(2)*(N0/2 + Sshot(Prx)/2 + Srin(RIN, Prx)/2) + 1/Nc*2*varQ));

%% Figures
Colors = num2cell(lines(6), 2);
figure, hold on, grid on, box on
plot(-10, -10, '--k', 'LineWidth', 2)
plot(-10, -10, '-k', 'LineWidth', 2)
plot(PrxdBm, SNR16QAM, '--b', 'LineWidth', 2)
plot(PrxdBm, SNR16QAMQ, '-b', 'LineWidth', 2)
plot(PrxdBm, SNR64QAM, '--r', 'LineWidth', 2)
plot(PrxdBm, SNR64QAMQ, '-r', 'LineWidth', 2)
legend('Without Quantization', 'With Quantization', 'Location', 'NorthWest')
xlabel('Received Power (dBm)', 'FontSize', 12)
ylabel('SNR (dB)',  'FontSize', 12)
axis([-16 8 0 40])
set(gca, 'xtick', -16:2:8)

matlab2tikz('tikz\SNR_vs_Prec.tikz')



