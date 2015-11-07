clear, clc, close all

addpath data/
addpath ../../mpam
addpath ../../f
addpath ../f

M = 4;
ka = [0.1, 0.25];
BW = [20 100; 20 300]; % (10:2.5:50)*1e9;

ff = figure;  
legs = {};
count = 1;
for n = 1:length(ka)
    for m = 1:size(BW, 1)
        BW0GHz = BW(m, 1);
        GainBWGHz = BW(m, 2);
        S = load(sprintf('margin_vs_gain_results_%d-PAM_%s_ka=%d_BW0=%d_GainBW=%d', M, 'equally-spaced', round(100*ka(n)), BW0GHz, GainBWGHz));
        
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
        
        figure(ff), hold on, box on
        if m <= 3
            hlines(count) = plot(S.Gains, S.MargindB, '-');
            plot(S.Gopt_margin, S.OptMargindB, 'o', 'Color', get(hlines(m), 'Color'));
        else
            hlines(count) = plot(S.Gains, S.MargindB, '--');
            plot(S.Gopt_margin, S.OptMargindB, 's', 'Color', get(hlines(m), 'Color'));
        end
        
        legs = [legs sprintf('ka = %.2f, BW0 = %.2f, GainBW = %.2f', ka(n), BW0GHz, GainBWGHz)];
        count = count + 1;
        
        %% Get data
        Gopt(n, m) = S.Gopt_margin;
        S.apdG.Gain = Gopt(n, m);
        BWopt(n, m) = S.apdG.BW/1e9;
        Marginopt(n, m) = S.OptMargindB;
        
        G1dB(n, m) = interp1(S.MargindB(S.Gains <= S.Gopt_margin), S.Gains(S.Gains <= S.Gopt_margin), S.OptMargindB-1);
        S.apdG.Gain = G1dB(n, m);
        BW1dB(n, m) = S.apdG.BW/1e9;
        
        G2dB(n, m) = interp1(S.MargindB(S.Gains <= S.Gopt_margin), S.Gains(S.Gains <= S.Gopt_margin), S.OptMargindB-2);
        S.apdG.Gain = G2dB(n, m);
        BW2dB(n, m) = S.apdG.BW/1e9;
        
        fprintf('ka = %.2f, BW0 = %.2f, GainBW = %.2f\n', ka(n), BW0GHz, GainBWGHz)
        [Gopt(n, m), BWopt(n, m), Marginopt(n, m); G1dB(n, m) BW1dB(n, m), Marginopt(n, m)-1]
        
    end
end

figure(ff)
xlabel('APD Gain (Linear Units)')
ylabel('Margin Improvement (dB)')
leg = legend(hlines, legs);
set(leg, 'Location', 'SouthEast')
drawnow
