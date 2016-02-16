clear, clc, close all

addpath data/
addpath ../../mpam
addpath ../../f
addpath ../f
addpath ../
addpath ../../other

m2tikz = matlab2tikz();

M = 8;
ka = 0.1;
BW = [20 300]; % (10:2.5:50)*1e9;
lineStyle = {'-', '--'};
level_spacing = {'equally-spaced', 'optimized'};
modBWGHz = 30;
Gains = 1:0.5:20;
ReferencePowerdBm = -12.890616516961265; 
% power required to achieve 1.8e-4 with 4-PAM in an ideal channel with 
% input referred noise of 30 pA/sqrt(Hz)
Colors = {'blue', 'red', 'green'};

ff = figure;  
legs = {};
count = 1;
for n = 1:length(ka)
    for kk = 1:length(level_spacing)
        for m = 1:size(BW, 1)
            BW0GHz = BW(m, 1);
            GainBWGHz = BW(m, 2);
            S = load(sprintf('data/margin_vs_gain_%d-PAM_%s_ka=%d_BW0=%d_GainBW=%d_modBW=%d',...
                M, level_spacing{kk}, round(100*ka(n)), BW0GHz, GainBWGHz, modBWGHz));

            %% BER plot
            figure, hold on, box on
            leg = {};
            for k = 1:length(S.Gains)
                hline(k) = plot(S.tx.PtxdBm, log10(S.BER(k).gauss), '-');
                plot(S.tx.PtxdBm, log10(S.BER(k).awgn), '--', 'Color', get(hline(k), 'Color'))
                plot(S.tx.PtxdBm, log10(S.BER(k).count), '-o', 'Color', get(hline(k), 'Color'))
                leg = [leg sprintf('Gain = %.2f', S.Gains(k))];
            end
            xlabel('Received Power (dBm)')
            ylabel('log(BER)') 
            legend(hline, leg);
            axis([S.tx.PtxdBm(1) S.tx.PtxdBm(end) -8 0])
            set(gca, 'xtick', S.tx.PtxdBm)
            title(sprintf('%d-PAM, %s, ka=%.2f, BW0=%d, GainBW, %d, modBW=%d', M, level_spacing{kk}, ka(n), BW0GHz, GainBWGHz, modBWGHz));     

            %% Sensitivity vs gain plot
            % SensitivitydB = interp1(S.Gains, S.SensitivitydB, Gains, 'spline');
            
            SensitivitydB = ReferencePowerdBm - S.PrxdBm_BERtarget;
            OptSensitivitydB = ReferencePowerdBm - S.PrxdBm_BERtarget_opt;
            
            ind = (S.Gains == GainBWGHz/BW0GHz | S.Gains == GainBWGHz/BW0GHz-0.5 | S.Gains == GainBWGHz/BW0GHz+0.5);
            
            SensitivitydB(ind) = [];
            S.Gains(ind) = [];
            SensitivitydB = interp1(S.Gains, SensitivitydB, Gains, 'spline');
            
            if abs(max(SensitivitydB)-OptSensitivitydB) > 0.1
                warning('Incorret optimal gain found in file: Sensitivity_vs_gain_%d-PAM_%s_ka=%d_BW0=%d_GainBW=%d_modBW=%d',...
                M, level_spacing{kk}, round(100*ka(n)), BW0GHz, GainBWGHz, modBWGHz)
                disp('Expected: ')
                disp(max(SensitivitydB))
                disp('Found:')
                disp(OptSensitivitydB)
            end
            
            figure(ff), hold on, box on            
            hlines(count) = plot(Gains, SensitivitydB, lineStyle{kk});
            plot(S.Gopt_margin, OptSensitivitydB, 'o', 'Color', get(hlines(count), 'Color'));
            
            if strcmp(level_spacing{kk}, 'equally-spaced')
                label = sprintf('$k_A = %.1f$', ka(n));
                mark = 'o';
            else
                label = '';
                mark = 'o';
            end
            m2tikz.addplot(Gains, SensitivitydB, lineStyle{kk}, Colors{n}, 'none', label)
            m2tikz.addplot(S.Gopt_margin, OptSensitivitydB, 'none', Colors{n}, mark)                

            legs = [legs sprintf('ka = %.1f', ka(n))];
            count = count + 1;
        end
    end
end

figure(ff)
xlabel('APD Gain (Linear Units)')
ylabel('Sensitivity Improvement (dB)')
leg = legend(hlines, legs);
set(leg, 'Location', 'SouthEast')
drawnow

m2tikz.extract(gca, 'just axis'); % extract only axis, xlabel, etc
m2tikz.write('ISI_8PAM.tex')
