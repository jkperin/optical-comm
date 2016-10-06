clear, clc, close all

addpath ../
addpath ../f/
addpath ../../f
addpath ../../mpam
addpath ../../apd
addpath ../../soa

BERtarget = 1.8e-4;
Delay = 250;
Npols = [1 2];
LineWidth = 50:50:1000;
totalLinewidth = 2*1e3*LineWidth;

Qpsk = QAM(4, 2*112e9);
sim.ModFormat = Qpsk;
sim.BERtarget = BERtarget;
csi = sqrt(2)/2;
R = 1;
q = 1.60217662e-19;

SNRdB2PrxdBm = @(SNRdB) 10*log10(10^(SNRdB/10)*2*q*Qpsk.Rs/(R*1e-3));
SNRdBref = SNRreq(BERtarget, Qpsk.M, Qpsk.type);
PrefdBm = SNRdB2PrxdBm(SNRdBref);

CPRmethod = {'logic', 'costas'};
Colors = {'b', 'r', 'g', 'm', 'y'};
Markers = {'--o', '-s', ':'};
PrxdBm = zeros(2, 2, length(LineWidth));
for CPR = 1:2
    for k = 1:2 %1:length(Npols)
        for kk = 1:length(LineWidth)
            try 
                S = load(sprintf('QPSK_Analog_BER_L=0km_lamb=1310nm_ModBW=30GHz_OPLL-%s_Npol=%d_linewidth=%dkHz_delay=%dps.mat', CPRmethod{CPR}, Npols(k), LineWidth(kk), Delay));
            catch e
                disp(e.message);
                continue;
            end
            
            BERcount = 0;
            counter = 0;
            for n = 1:S.sim.Realizations
                if all(S.BER(n).count < 0.1);
                    BERcount = BERcount + S.BER(n).count;
                    counter = counter + 1;
                end
            end
            
            BERcount = BERcount/counter;
            
            figure(1), clf, hold on, box on
            hline = plot(S.Tx.PlaunchdBm, log10(BERcount), '-o', 'DisplayName', sprintf('CPR = %s | Npols=%d | linewidth = %d', CPRmethod{CPR}, Npols(k), LineWidth(kk)));
%             plot(S.Tx.PlaunchdBm, log10(S.BER.theory), '-', 'DisplayName', sprintf('CPR = %s | Npols=%d | linewidth = %d', CPRmethod{CPR}, Npols(k), LineWidth(kk)), 'Color', get(hline, 'Color'))

            BERcount = log10(BERcount);
            idx = find(BERcount <= -2 & BERcount >= -5.5);
            f = fit(S.Tx.PlaunchdBm(idx).', BERcount(idx).', 'poly2');
            [PrxdBm(CPR, k, kk), ~, exitflag] = fzero(@(x) f(x) - log10(BERtarget), -28);
            hline = plot(S.Tx.PlaunchdBm, f(S.Tx.PlaunchdBm), '-', 'Color', get(hline, 'Color'));
            axis([S.Tx.PlaunchdBm([1 end]) -8 0])
            if exitflag ~= 1
                disp('Interpolation failed')
                exitflag
                PrxdBm(CPR, k, kk) = NaN;
            end
%             pause(0.5)
            drawnow
        end
        figure(100), hold on, box on
        plot(totalLinewidth/1e3, squeeze(PrxdBm(CPR, k, :))-PrefdBm, Markers{k}, 'Color', Colors{CPR},...
            'DisplayName', sprintf('CPR = %s | Npols=%d', CPRmethod{CPR}, Npols(k)))
    end
end

Delay = Delay + 8;
for Ncpr = 1:2
    for k = 1:length(totalLinewidth)
        [wnOpt(Ncpr, k), ~, SNRdBstart] = optimizePLL(csi, Delay*1e-12, totalLinewidth(k), Ncpr, sim, false);
        
        [SNRdBpen, ~, exitflag] = fzero(@(SNRdB)...
            log10(ber_qpsk_imperfect_cpr(SNRdB, phase_error_variance(csi, wnOpt(Ncpr, k), Ncpr, Delay*1e-12, totalLinewidth(k), SNRdB, Qpsk.Rs))) - log10(BERtarget), SNRdBstart);
        
        if exitflag ~= 1
            warning('SNRdB calculation finished with exitflag = %d', exitflag);
        end
        
        PrxreqdBmTheory(Ncpr, k) = SNRdB2PrxdBm(SNRdBpen);
    end
end

figure(100), hold on, box on
plot(totalLinewidth/1e3, PrxreqdBmTheory(1, :)-PrefdBm, 'k', 'LineWidth', 2, 'DisplayName', '1 pol (theory)')
plot(totalLinewidth/1e3, PrxreqdBmTheory(2, :)-PrefdBm, 'k', 'LineWidth', 2, 'DisplayName', '2 pol (theory)')
axis([100 2000 0 2])
legend('-DynamicLegend')
m2tikz = matlab2tikz(gca);
m2tikz.write('process_Npols_penalty.tikz')
