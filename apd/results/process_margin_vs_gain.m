clear, clc, close all

addpath data/
addpath ../../mpam
addpath ../../f
addpath ../f
addpath ../
addpath ../../other

m2tikz = matlab2tikz();

M = 4;
ka = [0.1, 0.2];
BW = [20 100; 20 300; 30 100; 30 300]; % (10:2.5:50)*1e9;
lineStyle = {'-', '-', '--', '--'};
level_spacing = {'equally-spaced'};
modBWGHz = 30;
Gains = 1:0.5:30;
att = 0.35
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

            %% Margin vs gain plot
            % MargindB = interp1(S.Gains, S.MargindB, Gains, 'spline');
            
            MargindB = ReferencePowerdBm - S.PtxdBm_BERtarget;
            OptMargindB = ReferencePowerdBm - S.PtxdBm_BERtarget_opt;
            
            ind = (S.Gains == GainBWGHz/BW0GHz | S.Gains == GainBWGHz/BW0GHz-0.5 | S.Gains == GainBWGHz/BW0GHz+0.5);
            
            MargindB(ind) = [];
            Gains = S.Gains;
            Gains(ind) = [];
            MargindB = interp1(Gains, MargindB, S.Gains, 'spline');
            
            if abs(max(MargindB)-OptMargindB) > 0.1
                warning('Incorret optimal gain found in file: margin_vs_gain_%d-PAM_%s_ka=%d_BW0=%d_GainBW=%d_modBW=%d',...
                M, level_spacing{kk}, round(100*ka(n)), BW0GHz, GainBWGHz, modBWGHz)
                disp('Expected: ')
                disp(max(MargindB))
                disp('Found:')
                disp(OptMargindB)
            end
            
            figure(ff), hold on, box on            
            hlines(count) = plot(S.Gains, MargindB, lineStyle{kk});
            plot(S.Gopt_margin, OptMargindB, 'o', 'Color', get(hlines(count), 'Color'));
            
            if strcmp(level_spacing{kk}, 'equally-spaced')
                label = sprintf('$k_A = %.1f$', ka(n));
                mark = 'o';
            else
                label = '';
                mark = 'o';
            end
            m2tikz.addplot(S.Gains, MargindB, lineStyle{kk}, Colors{n}, 'none', label)
            m2tikz.addplot(S.Gopt_margin, OptMargindB, 'none', Colors{n}, mark)                

            legs = [legs sprintf('ka = %.1f', ka(n))];
            count = count + 1;
        end
    end
end

figure(ff)
xlabel('APD Gain (Linear Units)')
ylabel('Margin Improvement (dB)')
leg = legend(hlines, legs);
set(leg, 'Location', 'SouthEast')
drawnow

m2tikz.extract(gca, 'just axis'); % extract only axis, xlabel, etc
m2tikz.write('ISI_4PAM_30GHz.tex')
