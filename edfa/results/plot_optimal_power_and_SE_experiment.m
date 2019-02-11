%% Figures displaying optimized power allocation and spectral efficiency
clear, clc, close all

addpath ../
addpath ../data
addpath ../f/
addpath ../../f/

folder = 'edf_6m_1dB_margin';
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

        lnm = S.Signal.lnm;
               

        %% Linear regime
        % Optimized power loading
        figure(1), hold on, box on
        plot(lnm, S.lin.S.PdBm)
        xlabel('Wavelength (nm)')
        ylabel('Optmized channel power (dBm)')
  
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

        figure(4), hold on, box on
        hplot = plot(lnm, nlin.num.SE);
        plot(lnm, nlin.approx.SE, '--', 'Color', get(hplot, 'Color'))
        xlabel('Wavelength (nm)')
        ylabel('Spectral efficiency (bit/s/Hz)')
        legend('Numerical', 'Approximated')
        
        SE(p) = sum(nlin.num.SE);
        SE_approx(p) = sum(nlin.approx.SE);
        
        drawnow
    catch e
        warning(e.message)
    end
end

figure, hold on, box on
plot(pumpPowermW, SE, '-o')
plot(pumpPowermW, SE_approx, '-x')
xlabel('Pump power (mW)')
ylabel('Spectral efficiency (bit/s/Hz)')
legend('Numerical', 'Approximated')