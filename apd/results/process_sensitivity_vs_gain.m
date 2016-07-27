clear, clc, close all

addpath data/
addpath ../../mpam
addpath ../../f
addpath ../f
addpath ../
addpath ../../other

files = ls('sensitivity_vs_gain\*.mat')
ReferencePowerdBm = -12.890616516961265; 

% M = 8;
% ka = 0.1;
% BW = [20 300]; % (10:2.5:50)*1e9;
% lineStyle = {'-', '--'};
% level_spacing = {'equally-spaced', 'optimized'};
% modBWGHz = 30;
% Gains = 1:0.5:20;

% % power required to achieve 1.8e-4 with 4-PAM in an ideal channel with 
% % input referred noise of 30 pA/sqrt(Hz)
% Colors = {'blue', 'red', 'green'};

Gains = 1:0.25:20;
ff = figure;  
legs = {};
count = 1;
for n = 1:size(files, 1)
    file = files(n, :);
    
    S = load(['sensitivity_vs_gain\' file]);
    
    figure, hold on, box on
    leg = {};
    for k = 1:length(S.Gains)
        hline(k) = plot(S.Tx.PtxdBm, log10(S.BER(k).enum), '-');
        plot(S.Tx.PtxdBm, log10(S.BER(k).awgn), '--', 'Color', get(hline(k), 'Color'))
        plot(S.Tx.PtxdBm, log10(S.BER(k).count), '-o', 'Color', get(hline(k), 'Color'))
        leg = [leg sprintf('Gain = %.2f', S.Gains(k))];
        
%         log10ber = log10(S.BER(k).enum);
%         PrxdBm_BERtarget(k) = interp1(log10ber(log10ber > -4.5), S.PrxdBm(log10ber > -4.5), log10(S.sim.BERtarget));
    end    
    xlabel('Received Power (dBm)')
    ylabel('log(BER)') 
    legend(hline, leg);
    axis([S.Tx.PtxdBm(1) S.Tx.PtxdBm(end) -8 0])
    set(gca, 'xtick', S.Tx.PtxdBm)
    title(file); 
   
    
%     [PrxdBm_BERtarget_opt, PrxdBm_BERtarget_pin, PrxdBm_BERtarget]
    
    SensitivitydB = ReferencePowerdBm - S.PrxdBm_BERtarget;
    OptSensitivitydB = ReferencePowerdBm - S.PrxdBm_BERtarget_opt;
    if S.mpam.optimize_level_spacing
        lineStyle = '--';
    else
        lineStyle = '-';
    end
             
    figure(ff), hold on, box on    
    hlines(count) = plot(Gains, spline(S.Gains, SensitivitydB, Gains), lineStyle);
    plot(S.Apd.Gain, OptSensitivitydB, 'o', 'Color', get(hlines(count), 'Color'));
    legs = [legs file];
    count = count + 1;  
    drawnow
end

figure(ff)
xlabel('APD gain (linear units)')
ylabel('Sensitivity Improvement (dB)')
leg = legend(hlines, legs);
set(leg, 'Location', 'SouthEast')
drawnow

m2tikz = matlab2tikz(gca); % extract only axis, xlabel, etc
m2tikz.write('ISI_4PAM_20GHz_Mod=30GHz.tex')
