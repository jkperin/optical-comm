%% 
clear, clc, close all

addpath ../
addpath ../f/
addpath ../../f/

folder = 'capacity_vs_pump_power_PdBm';
edf_type = 'principles_type3';
pumpWavelengthnm = 980;
pumpPowermW = [30:5:100 150:50:1000]; %%[30:5:100 150:50:250 275:25:400 450:50:1000];
Nspans = 286; % 317
spanLengthKm = 50; % 

clipped = {};
Pin = zeros(length(pumpWavelengthnm), length(pumpPowermW));
Pout = zeros(length(pumpWavelengthnm), length(pumpPowermW));
Lopt = zeros(length(pumpWavelengthnm), length(pumpPowermW));
SEnum = zeros(length(pumpWavelengthnm), length(pumpPowermW));
SEapprox = zeros(length(pumpWavelengthnm), length(pumpPowermW));
BW = zeros(length(pumpWavelengthnm), length(pumpPowermW));
PCE = zeros(length(pumpWavelengthnm), length(pumpPowermW));
PCEmax = zeros(length(pumpWavelengthnm), length(pumpPowermW));
maxP = zeros(length(pumpWavelengthnm), length(pumpPowermW));
for n = 1:length(pumpWavelengthnm)
    for p = 1:length(pumpPowermW)
        try
            filename = sprintf('%s/capacity_vs_pump_power_EDF=%s_pump=%dmW_%dnm_L=%d_x_%dkm.mat',...
                folder, edf_type, pumpPowermW(p), pumpWavelengthnm(n), Nspans, spanLengthKm);
            disp(filename)
            S = load(filename);
            
            real(S.SE)
            
            kopt = S.kopt;
            if any(S.Sopt{kopt}.P > 0.95*S.PonVec(kopt))
                fprintf('## Clipped Ppump = %d mW\n', pumpPowermW(p));
                a = [2 1];
                kopt2 = a(kopt);
                if abs(diff(real(S.SE))) < 25
                    disp('selected other option')
                    kopt = kopt2;
                else
                    fprintf('capacity difference = %2.f\n', abs(diff(real(S.SE))))
                    clipped = [clipped {pumpPowermW(p), S.PonVec(S.kopt)}];
                    continue
                end
            end
            
            figure(495), clf, hold on, box on
            Line = {'--b', '-k'};
            plot(S.Sopt{1}.wavelength*1e9, S.Sopt{1}.PdBm, Line{1+(kopt == 1)}, 'DisplayName', sprintf('%.2f | %.2f', sum(S.num{1}.SE), sum(S.approx{1}.SE)))
            plot(S.Sopt{1}.wavelength([1 end])*1e9, Watt2dBm(S.PonVec(1))*[1 1], ':k', 'DisplayName', 'Power cap 1')
            plot(S.Sopt{2}.wavelength*1e9, S.Sopt{2}.PdBm, Line{1+(kopt == 2)}, 'DisplayName', sprintf('%.2f | %.2f', sum(S.num{2}.SE), sum(S.approx{2}.SE)))
            plot(S.Sopt{1}.wavelength([1 end])*1e9, Watt2dBm(S.PonVec(2))*[1 1], ':k', 'DisplayName', 'Power cap 2')
            xlabel('Wavelength (nm)')
            ylabel('Power (dBm)')
            legend('-dynamiclegend', 'Location', 'SouthEast')
            title(sprintf('Pump power = %d mW', pumpPowermW(p)))
                       
            Eopt = S.Eopt{kopt};
            Sopt = S.Sopt{kopt};
            num_opt = S.num{kopt};
            approx_opt = S.approx{kopt};
            Pin(n, p) = sum(Sopt.P);
            Pout(n, p) = sum(dBm2Watt(Sopt.PdBm + num_opt.GaindB));
            Lopt(n, p) = Eopt.L;
            SEnum(n, p) = sum(num_opt.SE);
            SEapprox(n, p) = sum(approx_opt.SE);
            BW(n, p) = 0.1*S.df*sum(Sopt.P ~= 0)/12.5e9;
            onChs = (Sopt.P ~= 0);
            PCE(n, p) = (Pout(n, p)-Pin(n, p))/S.Pump.P;
            PCEmax(n, p) = S.Pump.wavelength./1550e-9;
            maxP(n, p) = max(Sopt.P);
                
            if S.exitflag{kopt} == 1
                disp('PSO converged')
            else
                disp('PSO did not converge')
            end
                
            figure(108+n)
            lamb = Sopt.wavelength*1e9;
            subplot(221), hold on, box on
            PdBm = Sopt.PdBm;
            PdBm(Sopt.PdBm < -40) = NaN;
            hplot = plot(lamb, PdBm);
            xlabel('Wavelength (nm)')
            ylabel('Input power (dBm)')
            xlim(lamb([1 end]))
            grid on

            subplot(222), hold on, box on
            plot(lamb, num_opt.GaindB, '-', 'Color', get(hplot, 'Color'), 'DisplayName', 'Numerical');
            plot(lamb, approx_opt.GaindB, '--', 'Color', get(hplot, 'Color'), 'DisplayName', 'Approximated')
            % plot(lamb, 10*log10(1 + S.Pump.P*S.Pump.wavelength./(Sopt.P.*Sopt.wavelength)), 'k') %% max gain
            xlabel('Wavelength (nm)')
            ylabel('Gain (dB)')
            xlim(lamb([1 end]))
            grid on

            subplot(223), hold on, box on
            plot(lamb, Watt2dBm(num_opt.Pase), '-', 'Color', get(hplot, 'Color'), 'DisplayName', 'Numerical');
            plot(lamb, Watt2dBm(approx_opt.Pase), '--', 'Color', get(hplot, 'Color'), 'DisplayName', 'Approximated')
            xlabel('Wavelength (nm)')
            ylabel('ASE (dBm)')
            xlim(lamb([1 end]))
            grid on

            subplot(224), hold on, box on
            num_opt.SE(num_opt.SE < 1e-3) = NaN;
            approx_opt.SE(approx_opt.SE < 1e-3) = NaN;
            plot(lamb, num_opt.SE, '-', 'Color', get(hplot, 'Color'), 'DisplayName', 'Numerical');
            plot(lamb, approx_opt.SE, '--', 'Color', get(hplot, 'Color'), 'DisplayName', 'Approximated')
            xlabel('Wavelength (nm)')
            ylabel('Spectral efficiency (bits/s/Hz)')
            axis([lamb([1 end]) 0 ceil(max([num_opt.SE approx_opt.SE]))])
            grid on
                
            figure(208), hold on, box on
            hplot = plot(lamb, num_opt.SNRdB);
            plot(lamb, approx_opt.SNRdB, '--', 'Color', get(hplot, 'Color'))
            xlabel('Wavelength (nm)')
            ylabel('SNR (dB)')
            axis([lamb([1 end]) 0 20])
            
            drawnow
            
        catch e
            disp(e.message)
            Lopt(n, p) = NaN;
            Pout(n, p) = NaN;
            SEnum(n, p) = NaN;
            SEapprox(n, p) = NaN;
        end
    end
    
    figure(1), hold on, box on
    plot(pumpPowermW, Lopt(n, :), 'LineWidth', 2, 'DisplayName', sprintf('%d nm', pumpWavelengthnm(n)), 'LineWidth', 2)
    xlabel('Pump power (mW)', 'FontSize', 12)
    ylabel('Optimal EDF length (m)', 'FontSize', 12)
    legend('-dynamiclegend')
    set(gca, 'FontSize', 12)
    
    figure(2), hold on, box on
    plot(pumpPowermW, 100*PCE(n, :), 'LineWidth', 2, 'DisplayName', sprintf('%d nm', pumpWavelengthnm(n)))
    plot(pumpPowermW, 100*PCEmax(n, :), 'k', 'LineWidth', 2, 'DisplayName', 'Theoretical limit')
    xlabel('Pump power (mW)', 'FontSize', 12)
    ylabel('Power conversion efficiency (%)', 'FontSize', 12)
    legend('-dynamiclegend', 'Location', 'SouthEast')
    set(gca, 'FontSize', 12)
       
    figure(3), hold on, box on
    Cnum = S.df*SEnum(n, :)/1e12;
    fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0,-Inf],...
               'Upper',[Inf,Inf, Inf],...
               'StartPoint',[1 1 0]);
    ft = fittype('a*log2(1 + b*x) + c','options',fo);
    Cfit = fit(pumpPowermW.', Cnum.', ft)
    hplt = plot(pumpPowermW, Cnum, 'LineWidth', 2, 'DisplayName', sprintf('%d nm (numerical)', pumpWavelengthnm(n)));
    plot(pumpPowermW, S.df*SEapprox(n, :)/1e12, '--', 'Color', get(hplt, 'Color'), 'LineWidth', 2, 'DisplayName', sprintf('%d nm (analytical)', pumpWavelengthnm(n)));
    Ppump = linspace(10, 1500);
    plot(Ppump, Cfit(Ppump), ':k', 'DisplayName', 'log fit')
    xlabel('Pump power (mW)', 'FontSize', 12)
    ylabel('Total capacity per fiber (Tb/s)', 'FontSize', 12)
%     axis([25 100 0 350])
    legend('-dynamiclegend', 'Location', 'SouthEast')
    set(gca, 'FontSize', 12)
    
    figure(4), box on
    plot(pumpPowermW, BW(n, :), 'LineWidth', 2, 'DisplayName', sprintf('%d nm', pumpWavelengthnm(n)));
    xlabel('Pump power (mW)', 'FontSize', 12)
    ylabel('Bandwidth (nm)', 'FontSize', 12)
    legend('-dynamiclegend', 'Location', 'SouthEast')
    set(gca, 'FontSize', 12)
    
    figure(5), box on, hold on
    plot(pumpPowermW, Watt2dBm(maxP(n, :)), 'LineWidth', 2);
    plot(pumpPowermW, Watt2dBm(1/70*980/1550*(1e-3*pumpPowermW/(8-1))), 'k', 'LineWidth', 2);
    xlabel('Pump power (mW)', 'FontSize', 12)
    ylabel('Maximum power (dBm)', 'FontSize', 12)
    legend('-dynamiclegend', 'Location', 'SouthEast')
    set(gca, 'FontSize', 12)
    
    
    figure(6), box on, hold on
    plot(pumpPowermW, Watt2dBm(Pout(n, :)), 'LineWidth', 2, 'DisplayName', 'Output power');
    plot(pumpPowermW, Watt2dBm(Pin(n, :)), 'LineWidth', 2, 'DisplayName', 'Input power');
    %plot(pumpPowermW, Watt2dBm(Pin(n, :) + pumpWavelengthnm(n)/1550*pumpPowermW*1e-3), 'k', 'LineWidth', 2,'DisplayName', 'Theoretical limit');
    xlabel('Pump power (mW)', 'FontSize', 12)
    ylabel('Power (dBm)', 'FontSize', 12)
    legend('-dynamiclegend', 'Location', 'SouthEast')
    set(gca, 'FontSize', 12)
    
    drawnow
end

V = 12e3;
I = V/(2*1*Nspans*spanLengthKm);
Ptot = V*I/2;
Pedfa = 0.9*Ptot/(2*Nspans);

figure, clf, hold on, box on
overhead = 1.5;
Ndim = 4:100;
for Ppump = [100 500 1000]
    for k = 1:length(Ndim)
        Cp(k) = Ndim(k)*Cfit(Ppump/(overhead*Ndim(k)));
    end
    plot(Ndim, Cp, 'Displayname', sprintf('Ppump = %d mW', Ppump), 'LineWidth',2)
    [Copt, idx] = max(Cp);
%     plot(Ndim(idx), Copt, '*r')
end

set(gca, 'FontSize', 12)
xlabel('Spatial dimensions', 'FontSize', 12)
ylabel('Capacity (Tb/s)', 'FontSize', 12)
legend('-dynamiclegend')
% xlim([30 1000])


figure, hold on, box on
hplt = plot(10*log10(pumpPowermW), Cnum, 'LineWidth', 2, 'DisplayName', sprintf('%d nm (numerical)', pumpWavelengthnm(n)));
plot(10*log10(pumpPowermW), S.df*SEapprox(n, :)/1e12, '--', 'Color', get(hplt, 'Color'), 'LineWidth', 2, 'DisplayName', sprintf('%d nm (analytical)', pumpWavelengthnm(n)));
xlabel('Pump power (dBm)', 'FontSize', 12)
ylabel('Total capacity per fiber (Tb/s)', 'FontSize', 12)
legend('-dynamiclegend', 'Location', 'SouthEast')
