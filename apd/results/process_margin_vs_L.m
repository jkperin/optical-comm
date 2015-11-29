clear, clc, close all

addpath data_fiber/
addpath ../../mpam
addpath ../../f
addpath ../f
addpath ../
addpath ../../other


m2tikz1 = matlab2tikz();
m2tikz2 = matlab2tikz();

M = 4;
ka = [0.1 0.2];
BW = [20 100; 20 300]; % (10:2.5:50)*1e9;
lineStyle = {'-', '--'};
level_spacing = {'equally-spaced', 'optimized'};
marker = {'o', 's'};
modBWGHz = 30;
Gains = 1:0.5:30;
pen = 1; % penalty of 1dB
Lkm = 0:10;
att = 0.35; % dB/km
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
                S = load(sprintf('data_fiber/Preq_vs_gain_L_%d-PAM_%s_ka=%d_BW0=%d_GainBW=%d_modBW=%d_L=%dkm',...
                M, level_spacing{kk}, round(100*ka(n)), BW0GHz, GainBWGHz, modBWGHz, Lkm(ll)));

                % BER plot
%                 figure, hold on, box on
%                 leg = {};
%                 for k = 1:length(S.Gains)
%                     hline(k) = plot(S.tx.PtxdBm, log10(S.BER(k).gauss), '-');
%                     plot(S.tx.PtxdBm, log10(S.BER(k).awgn), '--', 'Color', get(hline(k), 'Color'))
%                     plot(S.tx.PtxdBm, log10(S.BER(k).count), 'o', 'Color', get(hline(k), 'Color'))
%                     leg = [leg sprintf('Gain = %.2f', S.Gains(k))];
%                 end
%                 xlabel('Received Power (dBm)')
%                 ylabel('log(BER)') 
%                 legend(leg);
%                 axis([S.tx.PtxdBm(1) S.tx.PtxdBm(end) -8 0])
%                 set(gca, 'xtick', S.tx.PtxdBm)
%                 title(sprintf('%d-PAM, %s, ka = %.2f, BW0 = %.2f, GainBW', S.mpam.M, S.mpam.level_spacing, ka(n), BW0GHz, GainBWGHz))        
                
                %% Margin
                MargindB = ReferencePowerdBm - S.PrxdBm_BERtarget; % for all gains
                OptMargindB = ReferencePowerdBm - S.PrxdBm_BERtarget_opt; % for optimal gain

                ind = (S.Gains == GainBWGHz/BW0GHz | S.Gains == GainBWGHz/BW0GHz-0.5 | S.Gains == GainBWGHz/BW0GHz+0.5);

                MargindB(ind) = [];
                Gains = S.Gains;
                Gains(ind) = [];
                MargindB = interp1(Gains, MargindB, S.Gains, 'spline');
                
                MargindB_vs_L(ll, :) = MargindB;
                OptMargindB_vs_L(ll) = OptMargindB;
                Gopt_margin_vs_L(ll) = S.Gopt_margin;
                Gpen_margin_vs_L(ll) = interp1(MargindB, S.Gains, OptMargindB-pen);
            end
            
            figure(1000), hold on, box on
            plot(Lkm, OptMargindB_vs_L, lineStyle{kk}, 'Color', Colors{n}, 'Marker', marker{m})        
            legs = [legs sprintf('ka = %.2f, BW0 = %.2f GHz, GainBW=%.2f GHz', ka(n), BW0GHz, GainBWGHz)];
            m2tikz1.addplot(Lkm, OptMargindB_vs_L, lineStyle{kk}, Colors{n}, marker{m}, sprintf('ka = %.2f', ka(n)))
                        
            figure(1001), hold on, box on
            plot(Lkm, Gopt_margin_vs_L, lineStyle{kk}, 'Color', Colors{n}, 'Marker', marker{m}) 
            plot(Lkm, Gpen_margin_vs_L, lineStyle{kk}, 'Color', Colors{n}, 'Marker', marker{m}) 
            
            m2tikz2.addplot(Lkm, Gopt_margin_vs_L, lineStyle{kk}, Colors{n}, marker{m}, sprintf('ka = %.2f, GainBW=%.2f GHz'))
            m2tikz2.addplot(Lkm, Gpen_margin_vs_L, lineStyle{kk}, Colors{n}, marker{m}, sprintf('ka = %.2f, GainBW=%.2f GHz'))
            
            legs = [legs sprintf('ka = %.2f, BW0 = %.2f, GainBW', ka(n), BW0GHz, GainBWGHz)];
        end
    end
end
figure(1000)
xlabel('Fiber length (km)')
ylabel('Margin improvement (dB)')
axis([0 10 0 10])
legend(legs)
m2tikz1.extract(gca, 'just axis');
m2tikz1.write('fiber_4PAM_20GHz.tex');

figure(1001)
xlabel('Fiber length (km)')
ylabel('Optimal Gain (dB)')
legend(legs)

