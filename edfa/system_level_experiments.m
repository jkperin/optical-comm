%% System level experiments
clear, clc, close all

addpath results/
addpath f/
addpath ../f/

folder = 'results/';
filename = sprintf('%s/capacity_vs_pump_power_EDF=%s_pump=%dmW_%dnm_ChDf=%dGHz_L=%d_x_%dkm.mat',...
                folder, 'Corning (NEW)', 60, 980, 50, 220, 50);
disp(filename)
S = load(filename);

spanAttdB = S.spanAttdB;

OffPower = 1e-12;
df = S.problem.df;
E = S.nlin.E;
Pump = S.Pump;
Signal = S.nlin.S;
Signal.P(Signal.P == 0) = OffPower;
Signal_out = Signal;

OSNRinputdB = 45; 
ASEf = Channels(Signal.wavelength, dBm2Watt(Signal.PdBm - OSNRinputdB), 'forward');

PindBm = Signal.PdBm;
Namp = 5;
Nloop = ceil(S.Nspans/Namp);
loopLossdB = 6;

figure(1), hold on, box on
plot(Signal.lnm, Signal.PdBm)
xlabel('Wavelength (nm)')
ylabel('Input signal power (dBm)')
xlim([1520 1580])

for n = 1:Nloop
    fprintf('n = %d\n', n);
    ASEb = Channels(Signal.wavelength, 0, 'backward'); % zero backward ASE at every iteration
    
    for k = 1:Namp
        [GaindB, Ppump_out, Psignal_out, Pase, sol] = E.propagate(Pump, Signal, ASEf, ASEb, df, 'two-level', 50, false);
        
        Signal.P = dBm2Watt(Watt2dBm(Psignal_out) - spanAttdB);
        Signal.P(Signal.P < OffPower) = OffPower;
        ASEf.P = dBm2Watt(Watt2dBm(Pase) - spanAttdB);
        
        if any(sol.y(:) < 0)
            warning('EDF/propagate: solution contains negative power')
        end
    end
    
    [GaindB, Ppump_out, Psignal_out, Pase, sol] = E.propagate(Pump, Signal, ASEf, ASEb, df, 'two-level', 50, false);
    
    F = min(PindBm - Watt2dBm(Psignal_out) + loopLossdB, 0);
    
    Signal.P = dBm2Watt(Watt2dBm(Psignal_out) - loopLossdB + F );
    Signal.P(Signal.P < OffPower) = OffPower;
    ASEf.P = dBm2Watt(Watt2dBm(Pase) - loopLossdB + F);
                 
    figure(1), hold on, box on
    plot(Signal.lnm, Signal.PdBm)
    xlabel('Wavelength (nm)')
    ylabel('Input signal power (dBm)')
    xlim([1530 1565])
    
    figure(2), hold on, box on
    plot(Signal.lnm, GaindB)
    xlabel('Wavelength (nm)')
    ylabel('Gain (dB)')
    xlim([1530 1565])
    
    figure(3), hold on, box on
    plot(Signal.lnm, F)
    xlabel('Wavelength (nm)')
    ylabel('GFF (dB)')
    xlim([1530 1565])
                
    figure(5), hold on, box on
    plot(Signal.lnm, ASEf.PdBm)
    xlabel('Wavelength (nm)')
    ylabel('ASE power (dBm)')
    xlim([1530 1565])
        
    fprintf('Total input power = %.2f dBm\n', Watt2dBm(sum(Signal.P)))
    fprintf('Total ASE power = %.2f dBm\n', Watt2dBm(sum(ASEf.P)))
    drawnow
end

SNR = Signal.P./(ASEf.P + S.nlin.num.NL);
SE = 2*log2(1 + S.problem.Gap*SNR);

% figure, box on, hold on
% plot(Signal.lnm, 10*log10(SNR))
% plot(Signal.lnm, S.nlin_sfn.num.SNRdB)
% xlabel('Wavelength (nm)')
% ylabel('SNR (dB)')
% legend('Simulation', 'Predicted')
% xlim([1520 1580])

figure, box on, hold on
plot(Signal.lnm, SE)
plot(Signal.lnm, S.nlin.num.SE)
plot(Signal.lnm, S.nlin.approx.SE)
xlabel('Wavelength (nm)')
ylabel('Spectrale efficiency (bit/s/Hz)')
legend('Simulation', 'Predicted', 'Predicted approx')
% axis([1520 1580 0 4])
                        