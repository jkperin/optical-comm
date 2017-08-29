%% Validate classed EDF and channels
clear, clc, close all

% E = EDF(10, 'giles_ge:silicate');
E = EDF(4, 'principles_type3');
% E.plot('all');

Pump = Channels(1480e-9, 50e-3, 'forward');
Signal = Channels(linspace(1530, 1565, 50)*1e-9, 1e-4, 'forward');
ASEf = Channels(Signal.wavelength, 0, 'forward');
ASEb = Channels(Signal.wavelength, 0, 'backward');
% E.plot('coefficients', Signal.wavelength)

Lopt = E.optimal_length(Pump, Signal)

% figure, plot(Signal.wavelength, GaindB)
% figure, plot(Signal.wavelength, nsp)

% P = (1:20)*1e-3;
% for k = 1:length(P)
%     Pump.P =
% P(k);
%     GaindB(k) = 10*log10(E.analytical_gain(Pump, Signal));
%     nsp(k) = E.analytical_excess_noise(Pump, Signal);
% end

% figure, plot(P*1e3, GaindB)
% axis([P([1 end])*1e3 0 50])

% figure, plot(P*1e3, nsp)

dlamb = Signal.wavelength(2)-Signal.wavelength(1);
df = E.c/Signal.wavelength(1)-E.c/Signal.wavelength(2);
% df = 12.5e9;

GaindB_analytical = E.analytical_gain(Pump, Signal)
G_analytical = 10.^(GaindB_analytical/10);

Pase_analytical = E.analytical_ASE_PSD(Pump, Signal)*df;

[GaindB, Ppump_out, Psignal_out, Pase, sol] = E.two_level_system(Pump, Signal, ASEf, ASEb, df, 40)

figure, plot(Signal.wavelength*1e9, GaindB, Signal.wavelength*1e9, GaindB_analytical, '--')
figure, plot(Signal.wavelength*1e9, 10*log10(Pase/1e-3),...
    Signal.wavelength*1e9, 10*log10(Pase_analytical/1e-3), '--')
