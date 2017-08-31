%% Gain vs pump power
clear, clc, close all

Ppump = (0:10:80)*1e-3;

Pump = Channels(1480e-9, 50e-3, 'forward');
Signal = Channels(linspace(1530, 1565, 50)*1e-9, 1e-3, 'forward');
% Signal.P(1:15) = 0; % kills first and last 10 channels
% Signal.P(end:-1:end-15+1) = 0;
Signal.P(30:40) = 5e-3;
ASEf = Channels(Signal.wavelength, 0, 'forward');
ASEb = Channels(Signal.wavelength, 0, 'backward');

% E = EDF(8, 'giles_ge:silicate');
E = EDF(7, 'principles_type2');
Lopt = E.optimal_length(Pump, Signal)
% E.L = Lopt; 

dlamb = Signal.wavelength(2)-Signal.wavelength(1);
df = E.c/Signal.wavelength(1)-E.c/Signal.wavelength(2);

for k = 1:length(Ppump)
    Pump.P = Ppump(k);
    GaindB_semi_analytical(k, :) = E.semi_analytical_gain(Pump, Signal);
    [GaindB(k, :), Ppump_out, Psignal_out, Pase, sol] = E.two_level_system(Pump, Signal, ASEf, ASEb, df, 40);
    1;
end

% GaindB_anlaytical = E.analytical_gain(Pump, Signal);

% figure, box on, hold on
% for k = 1:5:50
%     hplot = plot(Ppump*1e3, GaindB(:, k), '-', 'DisplayName', sprintf('%.2f nm', Signal.wavelength(k)*1e9));
%     plot(Ppump*1e3, GaindB_semi_analytical(:, k), '--', 'Color', get(hplot, 'Color'), 'HandleVisibility','off');
% %     plot(Ppump([1 end])*1e3, GaindB_anlaytical(k)*[1 1], ':', 'Color', get(hplot, 'Color'), 'HandleVisibility','off');
% end
% legend('-DynamicLegend')
% xlabel('Pump power (mW)', 'FontSize', 12)
% ylabel('Gain (dB)', 'FontSize', 12)


figure(5), box on, hold on
for k = length(Ppump):-1:1
    hplot = plot(Signal.wavelength*1e9, GaindB(k, :), '-', 'DisplayName', sprintf('%1.f mW', Ppump(k)*1e3));
    plot(Signal.wavelength*1e9, GaindB_semi_analytical(k, :), '--', 'Color', get(hplot, 'Color'), 'HandleVisibility','off');
%     plot(Ppump([1 end])*1e3, GaindB_anlaytical(k)*[1 1], ':', 'Color', get(hplot, 'Color'), 'HandleVisibility','off');
end
legend('-DynamicLegend')
xlabel('Wavelength (nm)', 'FontSize', 12)
ylabel('Gain (dB)', 'FontSize', 12)
title(sprintf('EDF %s, L = %.1f m, %s pump', E.type, E.L, Pump.direction), 'Interpreter','none')

dGain = diff(GaindB, 1, 1);
lamb = repmat(Signal.wavelength*1e9, length(Ppump)-1, 1);
figure(5) 
quiver(lamb, GaindB(2:end, :), zeros(size(dGain)), dGain)
axis([1525 1570 -15 18])

    