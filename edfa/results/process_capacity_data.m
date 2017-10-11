%% 
clear, clc, close all

addpath ../
addpath ../f/

folder = 'capacity_vs_pump_power_EDF=principles_type3';
% _pump=5mW_980nm_L=280_x_50km'
edf_type = 'principles_type3';
pumpWavelengthnm = [980 1480];
pumpPowermW = 30:5:100;
Nspans = 280;
spanLengthKm = 50;

for n = 1:length(pumpWavelengthnm)
    for p = 1:length(pumpPowermW)
        try
            S = load(sprintf('%s/capacity_vs_pump_power_EDF=%s_pump=%dmW_%dnm_L=%d_x_%dkm.mat',...
                folder, edf_type, pumpPowermW(p), pumpWavelengthnm(n), Nspans, spanLengthKm));

            Lopt(n, p) = S.E.L;
            SEnum(n, p) = sum(S.num.SE);
            SEapprox(n, p) = sum(S.approx.SE);
            
            Signal = S.Signal;
            
            figure(108+n)
            subplot(221), hold on, box on
            plot(Signal.wavelength*1e9, Signal.PdBm)
            xlabel('Wavelength (nm)')
            ylabel('Power (dBm)')
            xlim(Signal.wavelength([1 end])*1e9)

            subplot(222), hold on, box on
            hplot = plot(Signal.wavelength*1e9, S.num.GaindB, 'DisplayName', 'Numerical');
            plot(Signal.wavelength*1e9, S.approx.GaindB, '--', 'Color', get(hplot, 'Color'), 'DisplayName', 'Approximated')
            xlabel('Wavelength (nm)')
            ylabel('Gain (dB)')
%             legend('-dynamicLegend', 'Location', 'Best')
            xlim(Signal.wavelength([1 end])*1e9)

            subplot(223), hold on, box on
            hplot = plot(Signal.wavelength*1e9, Watt2dBm(S.num.Pase), 'DisplayName', 'Numerical');
            plot(Signal.wavelength*1e9, Watt2dBm(S.approx.Pase), '--', 'Color', get(hplot, 'Color'), 'DisplayName', 'Approximated')
            xlabel('Wavelength (nm)')
            ylabel('ASE (dBm)')
%             legend('-dynamicLegend', 'Location', 'Best')
            xlim(Signal.wavelength([1 end])*1e9)

            subplot(224), hold on, box on
            S.num.SE(S.num.SE < 1e-3) = NaN;
            S.approx.SE(S.approx.SE < 1e-3) = NaN;
            hplot = plot(Signal.wavelength*1e9, S.num.SE, 'DisplayName', 'Numerical');
            plot(Signal.wavelength*1e9, S.approx.SE, '--', 'Color', get(hplot, 'Color'), 'DisplayName', 'Approximated')
            xlabel('Wavelength (nm)')
            ylabel('Spectral efficiency (bits/s/Hz)')
%             legend('-dynamicLegend', 'Location', 'Best')
            axis([Signal.wavelength([1 end])*1e9 0 ceil(max([S.num.SE S.approx.SE]))])
            drawnow
            
        catch e
            disp(e.message)
            Lopt(n, p) = NaN;
            SEnum(n, p) = NaN;
            SEapprox(n, p) = NaN;
        end
    end
    
    figure(1), hold on, box on
    plot(pumpPowermW, Lopt(n, :), 'DisplayName', sprintf('%d nm', pumpWavelengthnm(n)))
    xlabel('Pump power (mW)')
    ylabel('Optimal EDF length (m)')
    legend('-dynamiclegend')
       
    figure(2), hold on, box on
    hplt = plot(pumpPowermW, SEnum(n, :), 'DisplayName', sprintf('%d nm', pumpWavelengthnm(n)));
    plot(pumpPowermW, SEapprox(n, :), '--', 'Color', get(hplt, 'Color'), 'DisplayName', sprintf('%d nm', pumpWavelengthnm(n)));
    xlabel('Pump power (mW)')
    ylabel('Total spectrum efficiency (bits/s/Hz)')
    legend('-dynamiclegend')
    drawnow
end
