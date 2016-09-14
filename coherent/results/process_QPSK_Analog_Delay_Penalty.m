clear, clc, close all

BERtarget = 1e-3;
loopBandwidth = [50 100:100:1000];

filename = 'QPSK_Analog_Delay_Penalty_OPLL_costas_linewidth=200kHz_loopBandwidthMHz=';

for k = 1:length(loopBandwidth)
    S = load(sprintf('%s%d.mat', filename, loopBandwidth(k)));
    
    for kk = 1:length(S.Delayps)
        figure(k), hold on, box on
        h = plot(S.PlaunchdBm, log10(S.BER{kk}.count), '-o');
        plot(S.PlaunchdBm, log10(S.BERopt{kk}.count), '--o', 'Color', get(h, 'Color'))
        title(sprintf('loopBandwidth = %d', loopBandwidth(k)))
        
        BERcount = log10(S.BER{kk}.count);
        idx = not(isinf(BERcount));
        f = fit(S.PlaunchdBm(idx).', BERcount(idx).', 'linearinterp');
        [PrxdBm(k, kk), ~, exitflag] = fzero(@(x) f(x) - log10(BERtarget), -28);
        if exitflag ~= 1
            disp('Interpolation failed')
            exitflag
            PrxdBm(k, kk) = NaN;
        end
        
        BERcount = log10(S.BERopt{kk}.count);
        idx = not(isinf(BERcount));
        fopt = fit(S.PlaunchdBm(idx).', BERcount(idx).', 'linearinterp');
        [PrxdBmOpt(k, kk), ~, exitflag] = fzero(@(x) fopt(x) - log10(BERtarget), -28);
        
        if exitflag ~= 1
            disp('Interpolation failed')
            exitflag
            PrxdBmOpt(k, kk) = NaN;
        end
    end
    
    figure(100), hold on, box on
    hline = plot(S.Delayps, PrxdBm(k, :), '-o', 'DisplayName', num2str(loopBandwidth(k)));
    plot(S.Delayps, PrxdBmOpt(k, :), '--s', 'Color', get(hline, 'Color'));
    drawnow
    
end

figure(100)
legend('-DynamicLegend')
axis([0 400 -40 -20])



