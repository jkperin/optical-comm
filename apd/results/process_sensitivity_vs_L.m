clear, clc, close all

addpath ../../mpam
addpath ../../f
addpath ../f
addpath ../
addpath ../../other

m2tikz1 = matlab2tikz();
m2tikz2 = matlab2tikz();

M = 4;
ka = [0.1, 0.2];
BW = [20 100; 20 300]; % (10:2.5:50)*1e9;
lineStyle = {'-', '--'};
level_spacing = {'equally-spaced', 'optimized'};
marker = {'o', 's'};
modBWGHz = 30;
Lkm = 0:15;
ReferencePowerdBm = -12.890616516961265; 
% power required to achieve 1.8e-4 with 4-PAM in an ideal channel with 
% input referred noise of 30 pA/sqrt(Hz)
Colors = {'blue', 'red', 'green'};

legs = {};
for n = 1:length(ka)
    for kk = 1:length(level_spacing)
        for m = 1:size(BW, 1)
            for ll = 1:length(Lkm)
                BW0GHz = BW(m, 1);
                GainBWGHz = BW(m, 2);
                S = load(sprintf('fiber_1250nm/Preq_vs_gain_L_%d-PAM_%s_ka=%d_BW0=%d_GainBW=%d_modBW=%d_L=%dkm',...
                M, level_spacing{kk}, round(100*ka(n)), BW0GHz, GainBWGHz, modBWGHz, Lkm(ll)));

                % BER plot
                figure, hold on, box on
                leg = [];
                for k = 1:length(S.Gains)
                    hline(k) = plot(S.tx.PtxdBm, log10(S.BER(k).gauss), '-');
                    plot(S.tx.PtxdBm, log10(S.BER(k).awgn), '--', 'Color', get(hline(k), 'Color'))
                    plot(S.tx.PtxdBm, log10(S.BER(k).count), 'o', 'Color', get(hline(k), 'Color'))
                    leg = [leg sprintf('Gain = %.2f', S.Gains(k))];
                end
                hline(end+1) = plot(S.tx.PtxdBm, log10(S.ber_apd_opt.gauss), '-k');
                plot(S.tx.PtxdBm, log10(S.ber_apd_opt.awgn), '--', 'Color', get(hline(end), 'Color'))
                plot(S.tx.PtxdBm, log10(S.ber_apd_opt.count), 'o', 'Color', get(hline(end), 'Color'))
                leg = [leg sprintf('Optimal Gain = %.2f', S.Gopt)];
                xlabel('Received Power (dBm)')
                ylabel('log(BER)') 
                axis([S.tx.PtxdBm(1) S.tx.PtxdBm(end) -8 0])
                set(gca, 'xtick', S.tx.PtxdBm)
                title(sprintf('%d-PAM, %s, ka = %.2f, BW0 = %.2f, GainBW', S.mpam.M, S.mpam.level_spacing, ka(n), BW0GHz, GainBWGHz))        
                drawnow
                
                %% Margin
                SensitivityImprovementdB(ll) = ReferencePowerdBm - S.PrxdBm_BERtarget_opt; % for all gains 
                Gopt(ll) = S.Gopt;
            end
            
            figure(1000), hold on, box on
            plot(Lkm, SensitivityImprovementdB, lineStyle{kk}, 'Color', Colors{n}, 'Marker', marker{m})        
            legs = [legs sprintf('ka = %.2f, BW0 = %.2f GHz, GainBW=%.2f GHz', ka(n), BW0GHz, GainBWGHz)];
            m2tikz1.addplot(Lkm, SensitivityImprovementdB, lineStyle{kk}, Colors{n}, marker{m}, sprintf('ka = %.2f', ka(n)))
                        
            figure(1001), hold on, box on
            plot(Lkm, Gopt, lineStyle{kk}, 'Color', Colors{n}, 'Marker', marker{m})             
            m2tikz2.addplot(Lkm, Gopt, lineStyle{kk}, Colors{n}, marker{m}, sprintf('ka = %.2f, GainBW=%.2f GHz'))
            legs = [legs sprintf('ka = %.2f, BW0 = %.2f, GainBW', ka(n), BW0GHz, GainBWGHz)];
        end
    end
end
figure(1000)
xlabel('Fiber length (km)')
ylabel('Sensitivity improvement (dB)')
axis([0 15 0 10])
legend(legs)
m2tikz1.extract(gca, 'just axis');
m2tikz1.write('fiber_4PAM_20GHz.tex');

figure(1001)
xlabel('Fiber length (km)')
ylabel('Optimal Gain (dB)')
legend(legs)

