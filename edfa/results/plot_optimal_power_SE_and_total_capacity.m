%% Figures displaying optimized power allocation and spectral efficiency
clear, clc, close all

addpath ../
addpath ../data
addpath ../f/
addpath ../../f/

folder = 'edf_6m_1dB_margin';
folder2 = 'edf_length_optimized';
% folder = 'edf_length_optimized';
edf_type = 'corning high NA';
pumpWavelengthnm = 980;
pumpPowermW = [40 60 80]; 
Nspans = 220; 
ChDf = 50;
spanLengthKm = 50; % 

for p = 1:length(pumpPowermW)
    try
        filename = sprintf('%s/capacity_vs_pump_power_EDF=%s_pump=%dmW_%dnm_ChDf=%dGHz_L=%d_x_%dkm.mat',...
            folder, edf_type, pumpPowermW(p), pumpWavelengthnm, ChDf, Nspans, spanLengthKm);
        disp(filename)
        S = load(filename);
        
        filename2 = sprintf('%s/capacity_vs_pump_power_EDF=%s_pump=%dmW_%dnm_ChDf=%dGHz_L=%d_x_%dkm.mat',...
            folder2, edf_type, pumpPowermW(p), pumpWavelengthnm, ChDf, Nspans, spanLengthKm);
        disp(filename2)
        S_Lopt = load(filename2);

        lnm = S.Signal.lnm;
  
        % Optimized power loading
        figure(1), hold on, box on
        hplot = plot(lnm, S.lin.S.PdBm, 'LineWidth', 2, 'DisplayName', sprintf('Pump = %d mW', pumpPowermW(p)));
        xlabel('Wavelength (nm)', 'FontSize', 12)
        ylabel('Optmized channel power (dBm)', 'FontSize', 12)
        set(gca, 'FontSize', 12)
        legend('-dynamiclegend')

        % Spectral efficiency
        figure(2), hold on, box on
        plot(lnm, S.lin.num.SE, '-', 'LineWidth', 2, 'Color', get(hplot, 'Color'), 'DisplayName', sprintf('Pump = %d mW', pumpPowermW(p)));
        plot(lnm, S.lin.approx.SE, '--', 'LineWidth', 2, 'Color', get(hplot, 'Color'), 'HandleVisibility', 'off')
        xlabel('Wavelength (nm)', 'FontSize', 12)
        ylabel('Spectral efficiency (bits/s/Hz)', 'FontSize', 12)
        set(gca, 'FontSize', 12)
        legend('-dynamiclegend')
        
        SE_num(p) = sum(S.lin.num.SE);
        SE_approx(p) = sum(S.lin.approx.SE);
        
        SE_num_Lopt(p) = sum(S_Lopt.lin.num.SE);
        SE_approx_Lopt(p) = sum(S_Lopt.lin.approx.SE);
        Lopt(p) = S_Lopt.lin.E.L;
        
        drawnow
    catch e
        warning(e.message)
    end
end

% Total capacity
figure(3), hold on, box on
plot(pumpPowermW, 50 * SE_num / 1e3, '-o', 'LineWidth', 2)
plot(pumpPowermW, 50 * SE_num_Lopt / 1e3, '-s', 'LineWidth', 2)
xlabel('Pump power (mW)', 'FontSize', 12)
ylabel('Total capacity (Tb/s)', 'FontSize', 12)
legend('Chosen EDF length: L = 6 m', 'Optimal EDF length: L = 5.4 m', 'Location', 'SouthEast')
set(gca, 'FontSize', 12)
legend('-dynamiclegend')


Lopt

% figure(1)
% m = matlab2tikz(gca);
% m.write_tables('opt_power_linear_regime', 'same x')
% 
% figure(3)
% m = matlab2tikz(gca);
% m.write_tables('opt_SE_linear_regime', 'same x')
% 
% figure(2)
% m = matlab2tikz(gca);
% m.write_tables('opt_power_nonlinear_regime', 'same x')
% 
% figure(4)
% m = matlab2tikz(gca);
% m.write_tables('opt_SE_nonlinear_regime', 'same x')