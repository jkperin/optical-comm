%% Single channel test
clear, clc, close all

addpath data/
addpath f/
addpath ../f/

E = EDF(6, 'corning_type1');
% E.excess_loss =0.3;

Pump = Channels(980e-9, 60e-3, 'forward');
Signal = Channels(linspace(1535, 1565, 40)*1e-9, 0, 'forward');
ASEf = Channels(Signal.wavelength, 0, 'forward');
ASEb = Channels(Signal.wavelength, 0, 'backward');

PindBm = [-30:10];

for k = 1:length(PindBm)
    Signal.P = (dBm2Watt(PindBm(k))/40)*ones(size(Signal.P));
    [GdB, Ppump_out, Pout, Pase, sol] = E.propagate(Pump, Signal, ASEf, ASEb, 12.5e9, 'three-level', 50, false);
    
    GaindB(k, :) = GdB;
    Psignal_out(k, :) = Pout;
%     GaindB_sa(k) = E.semi_analytical_gain(Pump, Signal);
end

figure, plot(Watt2dBm(Psignal_out), GaindB)
ylabel('Gain (dB)')
xlabel('Output power (dBm)')

% figure, hold on, box on
% plot(PindBm, GaindB)
% plot(PindBm, GaindB_sa)

% figure, plot(PindBm, Watt2dBm(Psignal_out))