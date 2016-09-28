clear, clc, close all

addpath ../
addpath ../f/
addpath ../../f
addpath ../../mpam
addpath ../../apd
addpath ../../soa

BERtarget = 1.8e-4;
Npols = [1 2];
LineWidth = 0:100:1000;

CPRmethod = {'logic', 'costas'};
Colors = {'b', 'r', 'g', 'm', 'y'};
Markers = {'-', '--', ':'};
for CPR = 1:2
    for k = 1:2 %1:length(Npols)
        for kk = 1:length(LineWidth)
            counter = 1;
            try 
                S = load(sprintf('opll/QPSK_Analog_BER_OPLL_%s_Npol=%d_L=0km_linewidth=%dkHz_ideal=1_delay=200ps.mat', CPRmethod{CPR}, Npols(k), LineWidth(kk)));
            catch e
                disp(e.message);
                continue;
            end
            BERcount = S.BER.count;
            while true
                try
                    Snew = load(sprintf('opll/QPSK_Analog_BER_OPLL_%s_Npol=%d_L=0km_linewidth=%dkHz_ideal=1_delay=200ps(%d).mat', CPRmethod{CPR}, Npols(k), LineWidth(kk)), counter);
                    if Snew.BER.count < 0.1 % check if not cycle slip
                        BERcount = BERcount + Snew.BER.count;
                        counter = counter + 1;
                    end
                catch e
                    break;
                end
            end
            
            BERcount = BERcount/counter;
            
            figure(1), hold on, box on
            hline = plot(S.Tx.PlaunchdBm, log10(BERcount), '-o', 'DisplayName', sprintf('CPR = %s | Npols=%d | linewidth = %d', CPRmethod{CPR}, Npols(k), LineWidth(kk)));
%             plot(S.Tx.PlaunchdBm, log10(S.BER.theory), '-', 'DisplayName', sprintf('CPR = %s | Npols=%d | linewidth = %d', CPRmethod{CPR}, Npols(k), LineWidth(kk)), 'Color', get(hline, 'Color'))

            BERcount = log10(BERcount);
            idx = find(BERcount >= -4.5 & BERcount <= -2);
            f = fit(S.Tx.PlaunchdBm(idx).', BERcount(idx).', 'poly2');
            [PrxdBm(CPR, k, kk), ~, exitflag] = fzero(@(x) f(x) - log10(BERtarget), -28);
            hline = plot(S.Tx.PlaunchdBm, f(S.Tx.PlaunchdBm), '-', 'Color', get(hline, 'Color'));
            axis([S.Tx.PlaunchdBm([1 end]) -8 0])
            if exitflag ~= 1
                disp('Interpolation failed')
                exitflag
                PrxdBm(CPR, k, kk) = NaN;
            end
        end
        figure(100), hold on, box on
        plot(LineWidth, squeeze(PrxdBm(CPR, k, :)), Markers{CPR}, 'Color', Colors{k},...
            'DisplayName', sprintf('CPR = %s | Npols=%d', CPRmethod{CPR}, Npols(k)))
        drawnow
    end
end

axis([0 1000 -34 -28])
legend('-DynamicLegend')