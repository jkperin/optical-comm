%% 
clear, clc, close all

addpath ../
addpath ../f/
addpath ../../f/

folder = 'capacity_vs_pump_power_EDF=corning_type1_pump=980nm_L=286_x_50km';
edf_type = 'corning_type1';
pumpWavelengthnm = 980;
pumpPowermW = [30:5:150 200:100:500]; %%[30:5:100 150:50:250 275:25:400 450:50:1000];
Nspans = 286; % 317
spanLengthKm = 50; % 

clipped = {};
Pin = zeros(length(pumpWavelengthnm), length(pumpPowermW));
Pout = zeros(length(pumpWavelengthnm), length(pumpPowermW));
lin.Lopt = zeros(length(pumpWavelengthnm), length(pumpPowermW));
nlin.Lopt = zeros(length(pumpWavelengthnm), length(pumpPowermW));
nlin_local.Lopt = zeros(length(pumpWavelengthnm), length(pumpPowermW));
lin.SEnum = zeros(length(pumpWavelengthnm), length(pumpPowermW));
lin.SEapprox = zeros(length(pumpWavelengthnm), length(pumpPowermW));
nlin.SEnum = zeros(length(pumpWavelengthnm), length(pumpPowermW));
nlin.SEapprox = zeros(length(pumpWavelengthnm), length(pumpPowermW));
nlin_local.SEnum = zeros(length(pumpWavelengthnm), length(pumpPowermW));
nlin_local.SEapprox = zeros(length(pumpWavelengthnm), length(pumpPowermW));
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
            
            lnm = S.lamb*1e9;
            
            %% Linear regime
            % Optimized power loading
            figure(201), hold on, box on
            plot(lnm, S.lin.S.PdBm)
            xlabel('Wavelength (nm)')
            ylabel('Optmized channel power (dBm)')
            title('Linear regime')
            if S.lin.exitflag == 1
                disp('Linear regime PSO converged')
            else
                disp('Linear regime PSO did not converged')
            end
            
            figure(203), hold on, box on
            hplot = plot(lnm, S.lin.num.SE);
            plot(lnm, S.lin.approx.SE, '--', 'Color', get(hplot, 'Color'))
            xlabel('Wavelength (nm)')
            ylabel('Spectral efficiency (bit/s/Hz)')
            legend('Numerical', 'Approximated')
            title('Linear regime')
            
            lin.SEnum(n, p) = sum(S.lin.num.SE);
            lin.SEapprox(n, p) = sum(S.lin.approx.SE);
            lin.Lopt(n, p) = S.lin.E.L;          
            
            %% Nonlinear regime
            % Optimized power loading
            figure(202), hold on, box on
            hplot = plot(lnm, S.nlin.S.PdBm);
            plot(lnm, S.nlin_local.S.PdBm, '--', 'Color', get(hplot, 'Color'))
            legend('PSO', 'Hybrid optimization')
            xlabel('Wavelength (nm)')
            ylabel('Optmized channel power (dBm)')
            title('Nonlinear regime')
            ylim([-20 -10])
            
            if S.nlin.exitflag == 1
                disp('Nonlinear regime PSO converged')
            else
                disp('Nonlinear regime PSO did not converged')
            end
            
            if S.nlin_local.exitflag == 1
                disp('Nonlinear regime hybrid optimization converged')
            else
                disp('Nonlinear regime hybrid optimization did not converged')
            end
            
            figure(204), hold on, box on
            hplot = plot(lnm, S.nlin.num.SE);
            plot(lnm, S.nlin.approx.SE, '--', 'Color', get(hplot, 'Color'))
            hplot = plot(lnm, S.nlin_local.num.SE);
            plot(lnm, S.nlin_local.approx.SE, '--', 'Color', get(hplot, 'Color'))
            xlabel('Wavelength (nm)')
            ylabel('Spectral efficiency (bit/s/Hz)')
            legend('PSO: Numerical', 'PSO: Approximated', 'Hybrid: Numerical', 'Hybrid: Approximated')
            title('Nonlinear regime')
            
            nlin.SEnum(n, p) = sum(S.nlin.num.SE);
            nlin.SEapprox(n, p) = sum(S.nlin.approx.SE);
            nlin.Lopt(n, p) = S.nlin.E.L;
            
            nlin_local.SEnum(n, p) = sum(S.nlin_local.num.SE);
            nlin_local.SEapprox(n, p) = sum(S.nlin_local.approx.SE);
            nlin_local.Lopt(n, p) = S.nlin_local.E.L;
            
                            
%             subplot(222), hold on, box on
%             plot(lamb, num_opt.GaindB, '-', 'Color', get(hplot, 'Color'), 'DisplayName', 'Numerical');
%             plot(lamb, approx_opt.GaindB, '--', 'Color', get(hplot, 'Color'), 'DisplayName', 'Approximated')
%             % plot(lamb, 10*log10(1 + S.Pump.P*S.Pump.wavelength./(Sopt.P.*Sopt.wavelength)), 'k') %% max gain
%             xlabel('Wavelength (nm)')
%             ylabel('Gain (dB)')
%             xlim(lamb([1 end]))
%             grid on
% 
%             subplot(223), hold on, box on
%             plot(lamb, Watt2dBm(num_opt.Pase), '-', 'Color', get(hplot, 'Color'), 'DisplayName', 'Numerical');
%             plot(lamb, Watt2dBm(approx_opt.Pase), '--', 'Color', get(hplot, 'Color'), 'DisplayName', 'Approximated')
%             xlabel('Wavelength (nm)')
%             ylabel('ASE (dBm)')
%             xlim(lamb([1 end]))
%             grid on
% 
%             subplot(224), hold on, box on
%             num_opt.SE(num_opt.SE < 1e-3) = NaN;
%             approx_opt.SE(approx_opt.SE < 1e-3) = NaN;
%             plot(lamb, num_opt.SE, '-', 'Color', get(hplot, 'Color'), 'DisplayName', 'Numerical');
%             plot(lamb, approx_opt.SE, '--', 'Color', get(hplot, 'Color'), 'DisplayName', 'Approximated')
%             xlabel('Wavelength (nm)')
%             ylabel('Spectral efficiency (bits/s/Hz)')
%             axis([lamb([1 end]) 0 ceil(max([num_opt.SE approx_opt.SE]))])
%             grid on
%                 
%             figure(208), hold on, box on
%             hplot = plot(lamb, num_opt.SNRdB);
%             plot(lamb, approx_opt.SNRdB, '--', 'Color', get(hplot, 'Color'))
%             xlabel('Wavelength (nm)')
%             ylabel('SNR (dB)')
%             axis([lamb([1 end]) 0 20])
%             
            drawnow
            
        catch e
            warning(e.message)
            lin.Lopt(n, p) = NaN;
            lin.SEapprox(n, p) = NaN;
            lin.SEnum(n, p) = NaN;
            
            nlin.SEnum(n, p) = NaN;
            nlin.SEapprox(n, p) = NaN;
            nlin.Lopt(n, p) = NaN;
            
            nlin_local.SEnum(n, p) = NaN;
            nlin_local.SEapprox(n, p) = NaN;
            nlin_local.Lopt(n, p) = NaN;
        end
    end
    
    figure(1), hold on, box on
    plot(pumpPowermW, lin.Lopt(n, :), 'LineWidth', 2, 'DisplayName', 'Linear regime')
    plot(pumpPowermW, nlin.Lopt(n, :), 'LineWidth', 2, 'DisplayName', 'Nonlinear regime')
    plot(pumpPowermW, nlin_local.Lopt(n, :), 'LineWidth', 2, 'DisplayName', 'Nonlinear regime: hybrid opt')
    xlabel('Pump power (mW)', 'FontSize', 12)
    ylabel('Optimal EDF length (m)', 'FontSize', 12)
    legend('-dynamiclegend')
    set(gca, 'FontSize', 12)
    
%     figure(2), hold on, box on
%     plot(pumpPowermW, 100*PCE(n, :), 'LineWidth', 2, 'DisplayName', sprintf('%d nm', pumpWavelengthnm(n)))
%     plot(pumpPowermW, 100*PCEmax(n, :), 'k', 'LineWidth', 2, 'DisplayName', 'Theoretical limit')
%     xlabel('Pump power (mW)', 'FontSize', 12)
%     ylabel('Power conversion efficiency (%)', 'FontSize', 12)
%     legend('-dynamiclegend', 'Location', 'SouthEast')
%     set(gca, 'FontSize', 12)
       
    figure(3), hold on, box on
%     Cnum = S.df*SEnum(n, :)/1e12;
%     fo = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[0,0,-Inf],...
%                'Upper',[Inf,Inf, Inf],...
%                'StartPoint',[1 1 0]);
%     ft = fittype('a*log2(1 + b*x) + c','options',fo);
%     Cfit = fit(pumpPowermW.', Cnum.', ft)
%     hplt = plot(pumpPowermW, Cnum, 'LineWidth', 2, 'DisplayName', sprintf('%d nm (numerical)', pumpWavelengthnm(n)));
    hplot(1) = plot(pumpPowermW, S.df*lin.SEnum(n, :)/1e12, 'LineWidth', 2);
    hplot(2) = plot(pumpPowermW, S.df*nlin.SEnum(n, :)/1e12, 'LineWidth', 2);
    hplot(3) = plot(pumpPowermW, S.df*nlin_local.SEnum(n, :)/1e12, 'LineWidth', 2);
    plot(pumpPowermW, S.df*lin.SEnum(n, :)/1e12, '--', 'Color', get(hplot(1), 'Color'), 'LineWidth', 2);
    plot(pumpPowermW, S.df*nlin.SEnum(n, :)/1e12, '--', 'Color', get(hplot(2), 'Color'), 'LineWidth', 2);
    plot(pumpPowermW, S.df*nlin_local.SEnum(n, :)/1e12, '--', 'Color', get(hplot(3), 'Color'), 'LineWidth', 2);
    legend('Linear regime', 'Nonlinear regime', 'Nonlinear regime hybrid opt')
%     Ppump = linspace(10, 1500);
%     plot(Ppump, Cfit(Ppump), ':k', 'DisplayName', 'log fit')
    xlabel('Pump power (mW)', 'FontSize', 12)
    ylabel('Total capacity per fiber (Tb/s)', 'FontSize', 12)
%     axis([25 100 0 350])
%     legend('-dynamiclegend', 'Location', 'SouthEast')
    set(gca, 'FontSize', 12)
    
%     figure(4), box on
%     plot(pumpPowermW, BW(n, :), 'LineWidth', 2, 'DisplayName', sprintf('%d nm', pumpWavelengthnm(n)));
%     xlabel('Pump power (mW)', 'FontSize', 12)
%     ylabel('Bandwidth (nm)', 'FontSize', 12)
%     legend('-dynamiclegend', 'Location', 'SouthEast')
%     set(gca, 'FontSize', 12)
%     
%     figure(5), box on, hold on
%     plot(pumpPowermW, Watt2dBm(maxP(n, :)), 'LineWidth', 2);
%     plot(pumpPowermW, Watt2dBm(1/70*980/1550*(1e-3*pumpPowermW/(8-1))), 'k', 'LineWidth', 2);
%     xlabel('Pump power (mW)', 'FontSize', 12)
%     ylabel('Maximum power (dBm)', 'FontSize', 12)
%     legend('-dynamiclegend', 'Location', 'SouthEast')
%     set(gca, 'FontSize', 12)
%     
%     
%     figure(6), box on, hold on
%     plot(pumpPowermW, Watt2dBm(Pout(n, :)), 'LineWidth', 2, 'DisplayName', 'Output power');
%     plot(pumpPowermW, Watt2dBm(Pin(n, :)), 'LineWidth', 2, 'DisplayName', 'Input power');
%     %plot(pumpPowermW, Watt2dBm(Pin(n, :) + pumpWavelengthnm(n)/1550*pumpPowermW*1e-3), 'k', 'LineWidth', 2,'DisplayName', 'Theoretical limit');
%     xlabel('Pump power (mW)', 'FontSize', 12)
%     ylabel('Power (dBm)', 'FontSize', 12)
%     legend('-dynamiclegend', 'Location', 'SouthEast')
%     set(gca, 'FontSize', 12)
    
    drawnow
end

% V = 12e3;
% I = V/(2*1*Nspans*spanLengthKm);
% Ptot = V*I/2;
% Pedfa = 0.9*Ptot/(2*Nspans);
% 
% figure, clf, hold on, box on
% overhead = 1.5;
% Ndim = 4:100;
% for Ppump = [100 500 1000]
%     for k = 1:length(Ndim)
%         Cp(k) = Ndim(k)*Cfit(Ppump/(overhead*Ndim(k)));
%     end
%     plot(Ndim, Cp, 'Displayname', sprintf('Ppump = %d mW', Ppump), 'LineWidth',2)
%     [Copt, idx] = max(Cp);
% %     plot(Ndim(idx), Copt, '*r')
% end
% 
% set(gca, 'FontSize', 12)
% xlabel('Spatial dimensions', 'FontSize', 12)
% ylabel('Capacity (Tb/s)', 'FontSize', 12)
% legend('-dynamiclegend')
% % xlim([30 1000])
% 
% 
% figure, hold on, box on
% hplt = plot(10*log10(pumpPowermW), Cnum, 'LineWidth', 2, 'DisplayName', sprintf('%d nm (numerical)', pumpWavelengthnm(n)));
% plot(10*log10(pumpPowermW), S.df*SEapprox(n, :)/1e12, '--', 'Color', get(hplt, 'Color'), 'LineWidth', 2, 'DisplayName', sprintf('%d nm (analytical)', pumpWavelengthnm(n)));
% xlabel('Pump power (dBm)', 'FontSize', 12)
% ylabel('Total capacity per fiber (Tb/s)', 'FontSize', 12)
% legend('-dynamiclegend', 'Location', 'SouthEast')
