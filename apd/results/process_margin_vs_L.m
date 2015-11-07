clear, clc, close all

addpath data/
addpath ../../mpam
addpath ../../f
addpath ../f

M = 4;
ka = [0.1, 0.25];
BW = [20 100; 20 300]; % (10:2.5:50)*1e9;

ff = figure;  
ff2 = figure;
legs = {};
count = 1;
for n = 1:length(ka)
    for m = 1:size(BW, 1)
        BW0GHz = BW(m, 1);
        GainBWGHz = BW(m, 2);
        S = load(sprintf('margin_vs_L_results_%d-PAM_%s_ka=%d_BW0=%d_GainBW=%d', M, 'equally-spaced', round(100*ka(n)), BW0GHz, GainBWGHz));
        
        % BER plot
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
        legend(leg);
        axis([S.tx.PtxdBm(1) S.tx.PtxdBm(end) -8 0])
        set(gca, 'xtick', S.tx.PtxdBm)
        title(sprintf('%d-PAM, %s, ka = %.2f, BW0 = %.2f, GainBW', S.mpam.M, S.mpam.level_spacing, ka(n), BW0GHz, GainBWGHz))        
               
        %% Get data
        for k = 1:length(S.L)
            G1dB(k) = interp1(S.MargindB(k, S.Gains <= S.Gopt_margin(k)), S.Gains(S.Gains <= S.Gopt_margin(k)), S.OptMargindB(k)-1);
        end
        figure(ff), hold on, box on
        plot(S.L/1e3, S.OptMargindB, '-');
        legs = [legs sprintf('ka = %.2f, BW0 = %.2f, GainBW = %.2f', ka(n), BW0GHz, GainBWGHz)]; 
             
        figure(ff2), hold on, box on
        h = plot(S.L/1e3, S.Gopt_margin, '-');
        plot(S.L/1e3, G1dB, '--', 'Color', get(h, 'Color'));
    end
end

figure(ff)
xlabel('Fiber Length (km)')
ylabel('Margin Improvement (dB)')
legend(legs)

matlab2tikz('margin_vs_L_4PAM.tex')

figure(ff2)
xlabel('Fiber Length (km)')
ylabel('APD Gain (Linear Units)')
drawnow
