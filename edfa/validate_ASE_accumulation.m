%% Test ASE accumulation and power 
clear, clc, close all

addpath results/
addpath f/
addpath ../f/

folder = 'results/capacity_vs_pump_power_new';
filename = sprintf('%s/capacity_vs_pump_power_EDF=%s_pump=%dmW_%dnm_ChDf=%dGHz_L=%d_x_%dkm.mat',...
                folder, 'corning_type1', 30, 980, 50, 287, 50);
disp(filename)
S = load(filename);

spanAttdB = S.spanAttdB;

df = S.problem.df;
E = S.nlin_sfn.E;
Pump = S.Pump;
Signal = S.nlin_sfn.S;
Signal.P(Signal.P == 0) = eps;
ASEf = Channels(Signal.wavelength, 0, 'forward');
Signal_out = Signal;

for n = 1:S.Nspans
    fprintf('n = %d\n', n);
    ASEb = Channels(Signal.wavelength, 0, 'backward'); % zero backward ASE at every iteration
    
    figure(11), hold on, box on
    plot(Signal.lnm, Signal.PdBm)
    xlabel('Wavelength (nm)')
    ylabel('Signal power (dBm)')
    
    [GaindB, Ppump_out, Psignal_out, Pase, sol] = E.propagate(Pump, Signal, ASEf, ASEb, df, 'two-level', 50, false);
    
    F = min(-GaindB + spanAttdB, 0);
    
    if any(sol.y(:) < 0)
        warning('EDF/propagate: solution contains negative power')
    end
        
    figure(9), hold on, box on
    plot(Signal.lnm, GaindB)
    xlabel('Wavelength (nm)')
    ylabel('Gain (dB)')
    
    figure(10), hold on, box on
    plot(Signal.lnm, F)
    xlabel('Wavelength (nm)')
    ylabel('GFF (dB)')
        
    Signal.P = dBm2Watt(Watt2dBm(Psignal_out) + F - spanAttdB);
    Signal.P(Signal.P < eps) = eps;
    ASEf.P = dBm2Watt(Watt2dBm(Pase) + F - spanAttdB);
    
    figure(12), hold on, box on
    plot(Signal.lnm, ASEf.PdBm)
    xlabel('Wavelength (nm)')
    ylabel('ASE power (dBm)')
        
    fprintf('Total input power = %.2f dBm\n', Watt2dBm(sum(Signal.P)))
    fprintf('Total ASE power = %.2f dBm\n', Watt2dBm(sum(ASEf.P)))
    drawnow
end

SNR = Signal.P./(ASEf.P + S.nlin_sfn.num.NL);
SE = 2*log2(1 + S.problem.Gap*SNR);

figure, box on, hold on
plot(Signal.lnm, 10*log10(SNR))
plot(Signal.lnm, S.nlin_sfn.num.SNRdB)

figure, box on, hold on
plot(Signal.lnm, SE)
plot(Signal.lnm, S.nlin_sfn.num.SE)
            
Sin = S.nlin_sfn.S;
Sin.P = Sin.P + S.nlin_sfn.num.Pase;
ASEb = Channels(Signal.wavelength, 0, 'backward'); % zero backward ASE at every iteration
GaindB = E.propagate(Pump, Sin, ASEb, ASEb, df, 'two-level', 50, false);


            