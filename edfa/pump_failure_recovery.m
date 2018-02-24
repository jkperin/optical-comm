%% Check recovery from pump failure
clear, %close all

S = load('../results/Df=50GHz/capacity_vs_pump_power_EDF=corning_type1_pump=50mW_980nm_ChDf=50GHz_L=286_x_50km.mat');

Nit = 10;
k_fail = Inf; % amplifier with pump failure

Pump = S.Pump;
Pump_fail = Pump;
Pump_fail.P = 0*Pump_fail.P/2;

E = S.nlin_sfn.E;
Signal = S.nlin_sfn.S.sample(S.nlin_sfn.S.P ~= 0);
spanAttdB = S.problem.spanAttdB;
df = S.problem.df;

ASEf = Channels(Signal.wavelength, 0, 'forward');
ASEb = Channels(Signal.wavelength, 0, 'backward');

[GaindB, Ppump_out, SignalOut.P, Pase, sol] = E.two_level_system(Pump, Signal, ASEf, ASEb, df, 100, false);

GFFdB = spanAttdB - GaindB; % Gain flattening filter shape
GFF = 10.^(GFFdB/10);
assert(all(GFF < 1 & GFF > 0), 'Invalid GFF')

figure
plot(Signal.lnm, GFFdB);

%
diffdB = zeros(Nit, Signal.N);
Signalk = Signal;
figure(101), hold on, box on
plot(Signalk.lnm, Signalk.PdBm, 'k')
xlabel('Wavelength (nm)')
ylabel('Input power (dBm)')
Pase_total = 0;
for k = 1:Nit
    fprintf('k = %d\n', k);
    
    ASEf = Channels(Signal.wavelength, 0, 'forward');
    ASEb = Channels(Signal.wavelength, 0, 'backward');
    if k == k_fail
        [GaindB, Ppump_out, signalPout, Pase] = E.two_level_system(Pump_fail, Signalk, ASEf, ASEb, df, 100, false);
    else
        [GaindB, Ppump_out, signalPout, Pase] = E.two_level_system(Pump, Signalk, ASEf, ASEb, df, 100, false);
    end
    Pase_total = Pase_total + Pase;
    Signalk.P = signalPout + Pase;    
    Signalk.P = GFF.*10^(-spanAttdB/10).*Signalk.P;
    diffdB(k, :) = Signalk.PdBm - Signal.PdBm;
    
    figure(101), hold on, box on
    plot(Signalk.lnm, Signalk.PdBm)
    drawnow
end

figure(1)
plot((1:Nit)-k_fail, diffdB(:, [1 end]))
m = matlab2tikz(gca);
m.write_tables('recovery_from_pump_fail', 'same x')

figure(202), box on, hold on
plot(Signal.lnm, Watt2dBm(Pase_total))
