%% Normalize clipping and quantization noise
clear, clc, close all

r = linspace(2, 5);

%% Clip vs quantization

%% DC-OFDM (clipping twice symmetrically)
K = 1 - 2*qfunc(r);

% Normalized clipping noise (clipped symmetrically with the same clipping ratio)
Kr = @(r) (1 - K.^2 + 2*qfunc(r).*(r.^2 - 1) - 2*r/sqrt(2*pi).*exp(-r.^2/2))./K.^2;

% Normalized quantization noise (equivalent to quantization noise at the
% trnasmitter and receiver for DC-OFDM)
Kq = @(r, ENOB) K.*r.^2/(3*2^(2*ENOB))./K.^2;

figure1 = figure('Color',[1 1 1]);
hold on
plot(r, -10*log10(Kr(r)), ':k', 'LineWidth', 1.5)
plot(r, -10*log10(Kq(r, 5)), '--k', 'LineWidth', 1.5)
plot(r, -10*log10(Kr(r) + Kq(r, 5)), 'k', 'LineWidth', 1.5)
plot(r, -10*log10(Kq(r, 6)), '--k', 'LineWidth', 1.5)
plot(r, -10*log10(Kr(r) + Kq(r, 6)), 'k', 'LineWidth', 1.5)

legend('Clipping', 'Quant ENOB = 5', 'Quant ENOB = 6', 'FontSize', 14)
axis([r(1) r(end) 10 40])

legend('Clipping', 'Quantization', 'Clipping + Quantization');

% Create xlabel
xlabel('Clipping ratio','FontSize',14);

% Create ylabel
ylabel('SNR (dB)','FontSize',14);

set(gca, 'FontSize', 14)
set(gca, 'xtick', 2:0.5:5)
grid on
box on

% Create textbox
annotation(figure1,'textbox',...
    [0.74 0.4 0.2 0.1],...
    'String',{'ENOB = 6'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.74 0.6 0.2 0.1],...
    'String',{'ENOB = 5'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'LineStyle','none');

% saveas(gca, 'clipping_vs_quantization.emf')

% convert gca to latex
% matlab2tikz files must be in current directory
% matlab2tikz('tikz\clipping_vs_quantization.tikz')