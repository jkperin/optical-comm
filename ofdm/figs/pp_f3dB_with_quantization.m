%% Plot acoofdm of ACO-OFDM and DC-OFDM with quantization
clear, clc, close all

addpath data

%% Load data
load aco_ofdm_pp_quant
acoofdm = results;

load dc_ofdm_pp_quant
dcofdm = results;

%% CASE 1: ENOB = 5, Preemphasis, CS = 16
%% CASE 2: ENOB = 5, Palloc, CS = 16
%% CASE 3: ENOB = 6, Preemphasis, CS = 16
%% CASE 4: ENOB = 6, Preemphasis, CS = 64
%% CASE 5: ENOB = 6, Palloc, CS = 16
%% CASE 6:  ENOB = 6, Palloc, CS = 64

%% Figure 1: 16-QAM
% {ENOB = 5 (preemphasis), ENOB = 5 (palloc), ENOB = 6 (preemphasis), ENOB
% = 6 (palloc)}
markers = {'--ok', '--vk', '-ok', '-vk'};
Colors =  num2cell(lines(3), 2); 

figure1 = figure('Color',[1 1 1]);
hold on
plot(-10, -10, 'ok', -10, -10, 'vk', -10, -10, ':k', 'LineWidth', 2)

%% ACO-OFDM
% CASE 1: ENOB = 5, Preemphasis, CS = 16
plot(acoofdm(1).Fnl/1e9, acoofdm(1).pp_ook, markers{1}, 'Color', Colors{1}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1]);
% CASE 2: ENOB = 5, Palloc, CS = 16
plot(acoofdm(2).Fnl/1e9, acoofdm(2).pp_ook, markers{2}, 'Color', Colors{1}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
% CASE 3: ENOB = 6, Preemphasis, CS = 16
plot(acoofdm(3).Fnl/1e9, acoofdm(3).pp_ook, markers{3}, 'Color', Colors{1}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1]);
% CASE 5: ENOB = 6, Palloc, CS = 16
plot(acoofdm(5).Fnl/1e9, acoofdm(5).pp_ook, markers{4}, 'Color', Colors{1}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
% AWGN curve
plot([10 50], acoofdm(1).pp_awgn_ofdm*[1 1], ':', 'Color', Colors{1}, 'LineWidth', 2)

%% DC-OFDM
% CASE 1: ENOB = 5, Preemphasis, CS = 16
plot(dcofdm(1).Fnl/1e9, dcofdm(1).pp_ook, markers{1}, 'Color', Colors{3}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
% CASE 2: ENOB = 5, Palloc, CS = 16
plot(dcofdm(2).Fnl/1e9, dcofdm(2).pp_ook, markers{2}, 'Color', Colors{3}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
% CASE 3: ENOB = 6, Preemphasis, CS = 16
plot(dcofdm(3).Fnl/1e9, dcofdm(3).pp_ook, markers{3}, 'Color', Colors{3}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
% CASE 5: ENOB = 6, Palloc, CS = 16
plot(dcofdm(5).Fnl/1e9, dcofdm(5).pp_ook, markers{4}, 'Color', Colors{3}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
% AWGN curve
plot([10 50], dcofdm(1).pp_awgn_ofdm*[1 1], ':', 'Color', Colors{3}, 'LineWidth', 2)

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
    [0.70 0.23 0.165722821403953 0.0606893803546539],...
    'String',{'ACO-OFDM'},...
    'FontWeight','normal',...
    'FontSize',14,...
    'Color', Colors{1},...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.72 0.48 0.165722821403953 0.0606893803546539],...
    'String',{'DC-OFDM'},...
    'FontWeight','normal',...
    'FontSize',14,...
    'Color', Colors{3},...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.15 0.145 0.3 0.12],...
    'String',{'        ENOB = 5','        ENOB = 6'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'BackgroundColor',[1 1 1]);

% Create line
annotation(figure1,'line',[0.17 0.23],...
    [0.227 0.227],'LineStyle','--', 'LineWidth',1.5);

% Create line
annotation(figure1, 'line',[0.17 0.23],...
    [0.18 0.18],'LineWidth',1.5);


saveas(gca, 'pp_vs_f3dB_16QAM_quant', 'emf')

% convert gca to latex
% matlab2tikz files must be in current directory
matlab2tikz('tikz\pp_vs_f3dB_16QAM_quant.tikz')

%% Figure 1: 64-QAM
figure2 = figure('Color',[1 1 1]);
hold on
plot(-10, -10, 'ok', -10, -10, 'vk', -10, -10, ':k', 'LineWidth', 2)

%% ACO-OFDM
% CASE 4: ENOB = 6, Preemphasis, CS = 64
plot(acoofdm(4).Fnl/1e9, acoofdm(4).pp_ook, markers{3}, 'Color', Colors{1}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
% CASE 6:  ENOB = 6, Palloc, CS = 64
plot(acoofdm(6).Fnl/1e9, acoofdm(6).pp_ook, markers{4}, 'Color', Colors{1}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
% AWGN curve
plot([10 50], acoofdm(4).pp_awgn_ofdm*[1 1], ':', 'Color', Colors{1}, 'LineWidth', 2)

%% DC-OFDM
% CASE 4: ENOB = 6, Preemphasis, CS = 64
plot(dcofdm(4).Fnl/1e9, dcofdm(4).pp_ook, markers{3}, 'Color', Colors{3}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
% CASE 6:  ENOB = 6, Palloc, CS = 64
plot(dcofdm(6).Fnl/1e9, dcofdm(6).pp_ook, markers{4}, 'Color', Colors{3},'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [1 1 1])
% AWGN curve
plot([10 50], dcofdm(4).pp_awgn_ofdm*[1 1], ':', 'Color', Colors{3}, 'LineWidth', 2)

% 
legend('Const. bit loading & preemphasis', 'Opt. bit loading & power allocation', 'Ideal AWGN channel')

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
    [0.70 0.32 0.165722821403953 0.0606893803546539],...
    'String',{'ACO-OFDM'},...
    'FontWeight','normal',...
    'FontSize',14,...
    'Color', Colors{1},...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor',[1 1 1]);

% Create textbox
annotation(figure2,'textbox',...
    [0.7 0.66 0.165722821403953 0.0606893803546539],...
    'String',{'DC-OFDM'},...
    'FontWeight','normal',...
    'FontSize',14,...
    'Color', Colors{3},...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'BackgroundColor',[1 1 1]);


saveas(gca, 'pp_vs_f3dB_64QAM_quant', 'emf')

% convert gca to latex
% matlab2tikz files must be in current directory
matlab2tikz('tikz\pp_vs_f3dB_64QAM_quant.tikz')
