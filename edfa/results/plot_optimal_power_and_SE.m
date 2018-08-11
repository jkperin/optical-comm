%% Figures displaying optimized power allocation and spectral efficiency
clear, clc, close all

addpath ../
addpath ../data
addpath ../f/
addpath ../../f/

% folder = 'capacity_vs_pump_power_new';
% folder = 'capacity_vs_pump_NL_correct_NF_correct';
folder = 'capacity_vs_pump_power_high_na';
edf_type = 'corning high NA';
pumpWavelengthnm = 980;
pumpPowermW = [30 60 120 150]; 
Nspans = 287; 
ChDf = 50;
spanLengthKm = 50; % 

for p = 1:length(pumpPowermW)
    try
        filename = sprintf('%s/capacity_vs_pump_power_EDF=%s_pump=%dmW_%dnm_ChDf=%dGHz_L=%d_x_%dkm.mat',...
            folder, edf_type, pumpPowermW(p), pumpWavelengthnm, ChDf, Nspans, spanLengthKm);
        disp(filename)
        S = load(filename);

        lnm = S.Signal.lnm;
               

        %% Linear regime
        % Optimized power loading
        figure(1), hold on, box on
        plot(lnm, S.lin.S.PdBm)
        xlabel('Wavelength (nm)')
        ylabel('Optmized channel power (dBm)')
        ylim([-20 -10])

        figure(3), hold on, box on
        hplot = plot(lnm, S.lin.num.SE);
        plot(lnm, S.lin.approx.SE, '--', 'Color', get(hplot, 'Color'))
        xlabel('Wavelength (nm)')
        ylabel('Spectral efficiency (bit/s/Hz)')
        legend('Numerical', 'Approximated')   

        %% Nonlinear regime
%         nlin = S.nlin_unc;
        nlin = S.nlin_sfn;
        
        % Optimized power loading
        figure(2), hold on, box on
        plot(lnm, nlin.S.PdBm, '-')
        xlabel('Wavelength (nm)')
        ylabel('Optmized channel power (dBm)')
        ylim([-20 -10])

        figure(4), hold on, box on
        hplot = plot(lnm, nlin.num.SE);
        plot(lnm, nlin.approx.SE, '--', 'Color', get(hplot, 'Color'))
        xlabel('Wavelength (nm)')
        ylabel('Spectral efficiency (bit/s/Hz)')
        legend('Numerical', 'Approximated')
        
        
%         nlin.E.excess_loss = 0.2;
% %         nlin.E.L = nlin.E.L*1.08;
        [nsp, NFdB] = nlin.E.excess_noise(S.Pump, nlin.S, true);
  
        drawnow
    catch e
        warning(e.message)
    end
end

figure(1)
m = matlab2tikz(gca);
m.write_tables('opt_power_linear_regime', 'same x')

figure(3)
m = matlab2tikz(gca);
m.write_tables('opt_SE_linear_regime', 'same x')

figure(2)
m = matlab2tikz(gca);
m.write_tables('opt_power_nonlinear_regime', 'same x')

figure(4)
m = matlab2tikz(gca);
m.write_tables('opt_SE_nonlinear_regime', 'same x')