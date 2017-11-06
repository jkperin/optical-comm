%% 
clear, clc, close all

addpath ../
addpath ../f/
addpath ../../f/

folder = 'capacity_vs_pump_power_D=2';
% _pump=5mW_980nm_L=280_x_50km'
edf_type = 'principles_type3';
pumpWavelengthnm = [980 1480];
pumpPowermW = 30:10:100;
Nspans = 280;
spanLengthKm = 50;

for n = 1:length(pumpWavelengthnm)
    for p = 1:length(pumpPowermW)
        try
            filename = sprintf('%s/capacity_vs_pump_power_EDF=%s_pump=%dmW_%dnm_L=%d_x_%dkm.mat',...
                folder, edf_type, pumpPowermW(p), pumpWavelengthnm(n), Nspans, spanLengthKm);
            disp(filename)
            S = load(filename);

            Lopt(n, p) = S.E.L;
            S.num.SE = 2*S.num.SE; % 2 polarizations
            S.approx.SE = 2*S.approx.SE; 
            
            SEnum(n, p) = sum(S.num.SE);
            SEapprox(n, p) = sum(S.approx.SE);
            if S.exitflag == 1
                disp('PSO converged')
            else
                disp('PSO did not converge')
            end
                
%             [num, approx] = capacity_linear_regime(S.E, S.Pump, S.Signal, S.spanAttdB, S.Namp, S.df)
            
            for m = 1:3
                flat_fminbnd(m).num.SE(n, p) = sum(S.flat_fminbnd(m).num.SE);
                flat_fminbnd(m).approx.SE(n, p) = sum(S.flat_fminbnd(m).approx.SE);
                
                flat_interp(m).num.SE(n, p) = sum(S.flat_interp(m).num.SE);
                flat_interp(m).approx.SE(n, p) = sum(S.flat_interp(m).approx.SE);
            end
                        
            Signal = S.Signal;
            Ptot(n, p) = sum(Signal.P);
            
            figure(108+n)
            subplot(221), hold on, box on
            PdBm = Signal.PdBm;
            PdBm(Signal.PdBm < -40) = NaN;
            plot(Signal.wavelength*1e9, PdBm)
            xlabel('Wavelength (nm)')
            ylabel('Power (dBm)')
            xlim(Signal.wavelength([1 end])*1e9)
            grid on

            subplot(222), hold on, box on
            hplot = plot(Signal.wavelength*1e9, S.num.GaindB, 'DisplayName', 'Numerical');
            plot(Signal.wavelength*1e9, S.approx.GaindB, '--', 'Color', get(hplot, 'Color'), 'DisplayName', 'Approximated')
            xlabel('Wavelength (nm)')
            ylabel('Gain (dB)')
%             legend('-dynamicLegend', 'Location', 'Best')
            xlim(Signal.wavelength([1 end])*1e9)
            grid on

            subplot(223), hold on, box on
            hplot = plot(Signal.wavelength*1e9, Watt2dBm(S.num.Pase), 'DisplayName', 'Numerical');
            plot(Signal.wavelength*1e9, Watt2dBm(S.approx.Pase), '--', 'Color', get(hplot, 'Color'), 'DisplayName', 'Approximated')
            xlabel('Wavelength (nm)')
            ylabel('ASE (dBm)')
%             legend('-dynamicLegend', 'Location', 'Best')
            xlim(Signal.wavelength([1 end])*1e9)
            grid on

            subplot(224), hold on, box on
            S.num.SE(S.num.SE < 1e-3) = NaN;
            S.approx.SE(S.approx.SE < 1e-3) = NaN;
            hplot = plot(Signal.wavelength*1e9, S.num.SE, 'DisplayName', 'Numerical');
            plot(Signal.wavelength*1e9, S.approx.SE, '--', 'Color', get(hplot, 'Color'), 'DisplayName', 'Approximated')
            xlabel('Wavelength (nm)')
            ylabel('Spectral efficiency (bits/s/Hz)')
%             legend('-dynamicLegend', 'Location', 'Best')
            axis([Signal.wavelength([1 end])*1e9 0 ceil(max([S.num.SE S.approx.SE]))])
            grid on
                
            figure(208), hold on, box on
            hplot = plot(Signal.wavelength*1e9, S.num.SNRdB);
            plot(Signal.wavelength*1e9, S.approx.SNRdB, '--', 'Color', get(hplot, 'Color'))
            xlabel('Wavelength (nm)')
            ylabel('SNR (dB)')
            axis([Signal.wavelength([1 end])*1e9 0 15])
            
            drawnow
            
        catch e
            disp(e.message)
            Lopt(n, p) = NaN;
            SEnum(n, p) = NaN;
            SEapprox(n, p) = NaN;
        end
    end
    
    figure(1), hold on, box on
    plot(pumpPowermW, Lopt(n, :), 'DisplayName', sprintf('%d nm', pumpWavelengthnm(n)), 'LineWidth', 2)
    xlabel('Pump power (mW)', 'FontSize', 12)
    ylabel('Optimal EDF length (m)', 'FontSize', 12)
    legend('-dynamiclegend')
    set(gca, 'FontSize', 12)
    
    
    figure(2), hold on, box on
    plot(pumpPowermW, 100*(1e3*Ptot(n, :)./pumpPowermW), 'DisplayName', sprintf('%d nm', pumpWavelengthnm(n)))
    xlabel('Pump power (mW)', 'FontSize', 12)
    ylabel('Power conversion efficiency (%)', 'FontSize', 12)
    legend('-dynamiclegend')
       
    figure(3), hold on, box on
    hplt = plot(pumpPowermW, S.df*SEnum(n, :)/1e12, 'LineWidth', 2, 'DisplayName', sprintf('%d nm (numerical)', pumpWavelengthnm(n)));
    plot(pumpPowermW, S.df*SEapprox(n, :)/1e12, '--', 'Color', get(hplt, 'Color'), 'LineWidth', 2, 'DisplayName', sprintf('%d nm (analytical)', pumpWavelengthnm(n)));
%      plot(pumpPowermW, flat_fminbnd(1).num.SE(n, :), ':k', 'DisplayName', sprintf('%d nm', pumpWavelengthnm(n)));
%     plot(pumpPowermW, flat_fminbnd(2).num.SE(n, :), ':r', 'DisplayName', sprintf('%d nm', pumpWavelengthnm(n)));
    plot(pumpPowermW, 2*S.df*flat_fminbnd(3).num.SE(n, :)/1e12, ':', 'Color', get(hplt, 'Color'), 'LineWidth', 2, 'DisplayName', sprintf('%d nm flat power (numerical)', pumpWavelengthnm(n)));
%     plot(pumpPowermW, flat_interp(1).num.SE(n, :), ':k', 'DisplayName', sprintf('%d nm', pumpWavelengthnm(n)));
%     plot(pumpPowermW, flat_interp(2).num.SE(n, :), ':r', 'DisplayName', sprintf('%d nm', pumpWavelengthnm(n)));
%     plot(pumpPowermW, flat_interp(3).num.SE(n, :), ':', 'Color', get(hplt, 'Color'),  'DisplayName', sprintf('%d nm', pumpWavelengthnm(n)));
    xlabel('Pump power (mW)', 'FontSize', 12)
    ylabel('Total capacity per fiber (Tb/s)', 'FontSize', 12)
%     axis([25 100 0 350])
    legend('-dynamiclegend', 'Location', 'SouthEast')
    drawnow
end
