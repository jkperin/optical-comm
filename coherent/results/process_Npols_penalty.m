clear, clc, close all

addpath ../
addpath ../f/
addpath ../../f
addpath ../../mpam
addpath ../../apd
addpath ../../soa

BERtarget = 1e-3;
Npols = [1 2];
LineWidth = 0:200:2000;

CPRmethod = {'logic', 'costas'};

for CPR = 1:2
    for k = 1:length(Npols)
        for kk = 1:length(LineWidth)
%             S0 = load(sprintf('opll/QPSK_Analog_BER_OPLL_%s_Npol=%d_L=0km_linewidth=%dkHz_ideal=1_delay=0ps.mat', CPRmethod{CPR}, Npols(k), LineWidth(kk)));
            S1 = load(sprintf('opll/QPSK_Analog_BER_OPLL_%s_Npol=%d_L=0km_linewidth=%dkHz_ideal=1_delay=200ps.mat', CPRmethod{CPR}, Npols(k), LineWidth(kk)));

            
            S = S1;
%             BERcount = (31744*S0.BER.count(1:24)+64512*S1.BER.count)/(31744+64512);
            BERcount = S1.BER.count;
            
            figure(1), hold on, box on
            hline = plot(S.Tx.PlaunchdBm, log10(BERcount), '-o', 'DisplayName', sprintf('CPR = %s | Npols=%d | linewidth = %d', CPRmethod{CPR}, Npols(k), LineWidth(kk)));
%             plot(S.Tx.PlaunchdBm, log10(S.BER.theory), '-', 'DisplayName', sprintf('CPR = %s | Npols=%d | linewidth = %d', CPRmethod{CPR}, Npols(k), LineWidth(kk)), 'Color', get(hline, 'Color'))

            BERcount = log10(BERcount);
            idx = not(isinf(BERcount));
            f = fit(S.Tx.PlaunchdBm(idx).', BERcount(idx).', 'linearinterp');
            [PrxdBm(k, kk), ~, exitflag] = fzero(@(x) f(x) - log10(BERtarget), -28);
            if exitflag ~= 1
                disp('Interpolation failed')
                exitflag
                PrxdBm(k, kk) = NaN;
            end
        end
        figure(100), hold on, box on
        plot(LineWidth, PrxdBm(k, :), 'DisplayName', sprintf('CPR = %s | Npols=%d', CPRmethod{CPR}, Npols(k)))
    end
end

legend('-DynamicLegend')