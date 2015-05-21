clear, clc, close all

%% Power penalty vs cutoff frequency without quantization

data = {'pp_vs_f3dB_aco_ofdm_preemphasis_Nc64_CS16_Pb4';
        'pp_vs_f3dB_aco_ofdm_palloc_Nc64_CS16_Pb4';
        'pp_vs_f3dB_aco_ofdm_preemphasis_Nc64_CS64_Pb4';
        'pp_vs_f3dB_aco_ofdm_palloc_Nc64_CS64_Pb4'};

%% Figure 1
figure1 = figure;
hold on
for k = 1:length(data)
    load(data{k})
    
    eval([sprintf('plot%d = ', k) 'plot(' data{k} '.Fnl/1e9, ' data{k} '.power_pen_ook);'])
end

set(plot1, 'Marker', 'o', 'LineWidth', 2, 'Color', [0 0 0])
set(plot2, 'Marker', 'o', 'LineWidth', 2, 'Color', [1 0 0])
set(plot3, 'Marker', 's', 'LineWidth', 2, 'Color', [0 0 0])
set(plot4, 'Marker', 's', 'LineWidth', 2, 'Color', [1 0 0])

xlabel('Cut-off frequency (GHz)', 'FontSize', 12)
ylabel('Optical Power Penalty (dB)', 'FontSize', 12)
legend('Preemphasis & Constant Bit Loading', 'Optimal Power Alloc. & Bit Loading')
set(gca, 'FontSize', 12)

% Create textbox
% Create textbox
annotation(figure1,'textbox',...
    [0.753045009784736 0.160411312963517 0.138026418786693 0.0690476190476223],...
    'String',{'16-QAM'},...
    'FontWeight','bold',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.753201565557727 0.261485506990475 0.120012720156558 0.0690476190476224],...
    'String',{'64-QAM'},...
    'FontWeight','bold',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'LineStyle','none');


%% Figure 2
figure2 = figure;
hold on
for k = 1:length(data)
    load(data{k})
    
    eval([sprintf('plot%d = ', k) 'plot(' data{k} '.Fnl/1e9, ' data{k} '.PtxdBm);'])
end

set(plot1, 'Marker', 'o', 'LineWidth', 2, 'Color', [0 0 0])
set(plot2, 'Marker', 'o', 'LineWidth', 2, 'Color', [1 0 0])
set(plot3, 'Marker', 's', 'LineWidth', 2, 'Color', [0 0 0])
set(plot4, 'Marker', 's', 'LineWidth', 2, 'Color', [1 0 0])


xlabel('Cut-off frequency (GHz)', 'FontSize', 12)
ylabel('Transmitted power (dBn)', 'FontSize', 12)
legend('Preemphasis & Constant Bit Loading', 'Optimal Power Alloc. & Bit Loading')
set(gca, 'FontSize', 12)
title('NEP = 50 pA/Hz^{0.5}')

% Create textbox
% Create textbox
annotation(figure2,'textbox',...
    [0.753045009784736 0.160411312963517 0.138026418786693 0.0690476190476223],...
    'String',{'16-QAM'},...
    'FontWeight','bold',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure2,'textbox',...
    [0.753201565557727 0.261485506990475 0.120012720156558 0.0690476190476224],...
    'String',{'64-QAM'},...
    'FontWeight','bold',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'LineStyle','none');


%% Figure 3
figure3 = figure;
hold on
for k = [1 3]
    load(data{k})
    
    eval([sprintf('plot%d = ', k) 'plot(' data{k} '.Fnl/1e9, ' data{k} '.Fs);'])
end

set(plot1, 'Marker', 's', 'LineWidth', 2, 'Color', [0 0 0])
set(plot3, 'Marker', 'o', 'LineWidth', 2, 'Color', [1 0 0])

xlabel('Cut-off frequency (GHz)', 'FontSize', 12)
ylabel('Sampling rate (GHz)', 'FontSize', 12)
legend('16-QAM', '64-QAM')
set(gca, 'FontSize', 12)

