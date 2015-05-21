%% Plot acoofdm of ACO-OFDM and DC-OFDM without quantization
clear, clc, close all

addpath data

%% Load Data
load aco_ofdm_pp
acoofdm = results;

load dc_ofdm_pp
dcofdm = results;

load dc_ofdm_pp_full_dc
dcofdm_full = results;

%% 1. Preemphasis, CS = 16
%% 2. Preemphasis, CS = 64
%% 3. Palloc, CS = 16
%% 4. Palloc, CS = 64

markers = {'-o', '-v', ':'};
Colors = num2cell(lines(3), 2); %{'b', [0 .5 0], 'r'};

%% Figure 1: 16-QAM
figure1 = figure('Color',[1 1 1]);
hold on
plot(-10, -10, 'ok', -10, -10, 'vk', -10, -10, ':k', 'LineWidth', 2)

plot(acoofdm(1).Fnl/1e9, acoofdm(1).pp_ook, markers{1}, 'Color', Colors{1}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
plot(acoofdm(3).Fnl/1e9, acoofdm(3).pp_ook, markers{2}, 'Color', Colors{1}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
plot([10 50], acoofdm(1).pp_awgn_ofdm*[1 1], markers{3}, 'Color', Colors{1}, 'LineWidth', 2)

plot(dcofdm(1).Fnl/1e9, dcofdm(1).pp_ook, markers{1}, 'Color', Colors{2}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
plot(dcofdm(3).Fnl/1e9, dcofdm(3).pp_ook, markers{2}, 'Color', Colors{2}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
plot([10 50], dcofdm(1).pp_awgn_ofdm*[1 1], markers{3}, 'Color', Colors{2}, 'LineWidth', 2)

plot(dcofdm(1).Fnl/1e9, dcofdm_full(1).pp_ook, markers{1}, 'Color', Colors{3}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
plot(dcofdm(3).Fnl/1e9, dcofdm_full(3).pp_ook, markers{2}, 'Color', Colors{3}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])

%
leg1 = legend('Const. bit loading & preemphasis', 'Opt. bit loading & power allocation', 'Ideal AWGN channel');
set(leg1, 'Position',[0.232738095238095 0.727777777777778 0.666071428571429 0.187301587301587]);

% Create xlabel
xlabel('Cutoff frequency (GHz)','FontSize',14);

% Create ylabel
ylabel('Optical Power Penalty (dB)','FontSize',14);

set(gca, 'FontSize', 14)
axis([14.99 50 0 18])
grid on
box on

% Create textbox
annotation(figure1,'textbox',...
    [0.66 0.23 0.21 0.06],...
    'String',{'ACO-OFDM'},...
    'FontSize',14,...
    'Color', Colors{1},...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.66 0.415 0.24 0.04],...
    'String',{'r\sigma'' DC-OFDM'},...
    'FontSize',14,...
    'Color', Colors{2},...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.66 0.58 0.24 0.06],...
    'String',{'r\sigma DC-OFDM'},...
    'FontSize',14,...
    'Color', Colors{3},...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor',[1 1 1]);

saveas(gca, 'pp_vs_f3dB_16QAM_no_quant', 'png')

% convert gca to latex
% matlab2tikz files must be in current directory
matlab2tikz('tikz\pp_vs_f3dB_16QAM_no_quant.tikz')

%% Figure 2: 64-QAM
figure2 = figure('Color',[1 1 1]);
hold on
plot(-10, -10, 'ok', -10, -10, 'vk', -10, -10, ':k', 'LineWidth', 2)

plot(acoofdm(2).Fnl/1e9, acoofdm(2).pp_ook, markers{1}, 'Color', Colors{1}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
plot(acoofdm(4).Fnl/1e9, acoofdm(4).pp_ook, markers{2}, 'Color', Colors{1}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
plot([10 50], acoofdm(2).pp_awgn_ofdm*[1 1], ':', 'Color', Colors{1}, 'LineWidth', 1.5)

plot(dcofdm(2).Fnl/1e9, dcofdm(2).pp_ook, markers{1}, 'Color', Colors{2}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
plot(dcofdm(4).Fnl/1e9, dcofdm(4).pp_ook, markers{2}, 'Color', Colors{2}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
plot([10 50], dcofdm(2).pp_awgn_ofdm*[1 1], ':', 'Color', Colors{2}, 'LineWidth', 1.5)

plot(dcofdm(2).Fnl/1e9, dcofdm_full(2).pp_ook, markers{1}, 'Color', Colors{3}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
plot(dcofdm(4).Fnl/1e9, dcofdm_full(4).pp_ook, markers{2}, 'Color', Colors{3}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])

%
leg2 = legend('Const. bit loading & preemphasis', 'Opt. bit loading & power allocation', 'Ideal AWGN channel');
set(leg2, 'Position',[0.232738095238095 0.727777777777778 0.666071428571429 0.187301587301587]);

% Create xlabel
xlabel('Cutoff frequency (GHz)','FontSize',14);

% Create ylabel
ylabel('Optical Power Penalty (dB)','FontSize',14);

set(gca, 'FontSize', 14)
axis([14.99 50 0 18])
grid on
box on

% Create textbox
annotation(figure2,'textbox',...
    [0.67 0.28 0.22 0.06],...
    'String',{'ACO-OFDM'},...
    'FontSize',14,...
    'Color', Colors{1},...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor',[1 1 1]);

% Create textbox
annotation(figure2,'textbox',...
    [0.66 0.53 0.24 0.04],...
    'String',{'r\sigma'' DC-OFDM'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'Color', Colors{2},...
    'LineStyle','none',...
    'BackgroundColor',[1 1 1]);

% Create textbox
annotation(figure2,'textbox',[0.66 0.70 0.24 0.02],...
    'String',{'r\sigma DC-OFDM'},...
    'FontSize',14,...
    'Color', Colors{3},...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor',[1 1 1]);

saveas(gca, 'pp_vs_f3dB_64QAM_no_quant', 'png')

% convert gca to latex
% matlab2tikz files must be in current directory
matlab2tikz('tikz\pp_vs_f3dB_64QAM_no_quant.tikz')
