%% Compare numerical and analytical excess noise
clear, close all

addpath ../
addpath ../f/
addpath ../../f/

E = EDF(10, 'principles_type3');

df = 50e9;
dlamb = df2dlamb(df);
lamb = 1530e-9:dlamb:1565e-9;

Pon = 5e-5;
Signal = Channels(lamb, Pon, 'forward');
Pump = Channels(980e-9, 100e-3, 'forward');

Ppump = (50:50:200)*1e-3;
for k = 1:length(Ppump)
    Pump.P = Ppump(k);
    [nsp, NFdB] = E.excess_noise(Pump, Signal);
    nsp_analytical = E.analytical_excess_noise(Pump, Signal);
    mean(nsp)
    
    figure(1), hold on, box on
    hplot(k) = plot(Signal.wavelength*1e9, nsp, 'DisplayName', sprintf('%d mW', Ppump(k)*1e3));
    plot(Signal.wavelength*1e9, nsp_analytical, '--', 'Color', get(hplot(k), 'Color'))    
    drawnow
end

xlabel('Wavelength (nm)')
ylabel('Excess noise')
% legend('-dynamiclegend')
ylim([1 3])
