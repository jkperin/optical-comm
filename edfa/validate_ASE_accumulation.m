%% Test ASE accumulation and power 
clear, clc%, close all

addpath results/
addpath f/
addpath ../f/

folder = 'results/capacity_vs_pump_power_new';
filename = sprintf('%s/capacity_vs_pump_power_EDF=%s_pump=%dmW_%dnm_ChDf=%dGHz_L=%d_x_%dkm.mat',...
                folder, 'corning_type1', 30, 980, 50, 287, 50);
disp(filename)
S = load(filename);

spanAttdB = S.spanAttdB;

OffPower = 1e-12;
df = S.problem.df;
E = S.nlin_sfn.E;
Pump = S.Pump;
Signal = S.nlin_sfn.S;
Signal.P(Signal.P == 0) = OffPower;
ASEf = Channels(Signal.wavelength, 0, 'forward');
Signal_out = Signal;

for n = 1:S.Nspans
    fprintf('n = %d\n', n);
    ASEb = Channels(Signal.wavelength, 0, 'backward'); % zero backward ASE at every iteration
    
%     figure(1), hold on, box on
%     plot(Signal.lnm, Signal.PdBm)
%     xlabel('Wavelength (nm)')
%     ylabel('Input signal power (dBm)')
%     xlim([1520 1580])
    
    [GaindB, Ppump_out, Psignal_out, Pase, sol] = E.propagate(Pump, Signal, ASEf, ASEb, df, 'two-level', 50, false);
    
    F = min(-GaindB + spanAttdB, 0);
    
    if any(sol.y(:) < 0)
        warning('EDF/propagate: solution contains negative power')
    end
        
    figure(2), hold on, box on
    plot(Signal.lnm, GaindB)
    xlabel('Wavelength (nm)')
    ylabel('Gain (dB)')
    xlim([1520 1580])
    
    figure(3), hold on, box on
    plot(Signal.lnm, F)
    xlabel('Wavelength (nm)')
    ylabel('GFF (dB)')
    xlim([1520 1580])
    
    figure(4), hold on, box on
    plot(Signal.lnm, Watt2dBm(Psignal_out))
    xlabel('Wavelength (nm)')
    ylabel('Output signal power (dBm)')
    xlim([1520 1580])
        
    Signal.P = dBm2Watt(Watt2dBm(Psignal_out) + F - spanAttdB);
    Signal.P(Signal.P < OffPower) = OffPower;
    ASEf.P = dBm2Watt(Watt2dBm(Pase) + F - spanAttdB);
    
    figure(5), hold on, box on
    plot(Signal.lnm, ASEf.PdBm)
    xlabel('Wavelength (nm)')
    ylabel('ASE power (dBm)')
    xlim([1520 1580])
        
    fprintf('Total input power = %.2f dBm\n', Watt2dBm(sum(Signal.P)))
    fprintf('Total ASE power = %.2f dBm\n', Watt2dBm(sum(ASEf.P)))
    drawnow
end

figure(1)
m = matlab2tikz(gca);
m.write_tables('sim_input_power', 'same x')

figure(2)
m = matlab2tikz(gca);
m.write_tables('sim_gain', 'same x')

figure(3)
m = matlab2tikz(gca);
m.write_tables('sim_gff', 'same x')

figure(4)
m = matlab2tikz(gca);
m.write_tables('sim_output_power', 'same x')

figure(5)
m = matlab2tikz(gca);
m.write_tables('sim_ase', 'same x')

SNR = Signal.P./(ASEf.P + S.nlin_sfn.num.NL);
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
plot(Signal.lnm, S.nlin_sfn.num.SE)
plot(Signal.lnm, S.nlin_sfn.approx.SE)
xlabel('Wavelength (nm)')
ylabel('Spectrale efficiency (bit/s/Hz)')
legend('Simulation', 'Predicted', 'Predicted approx')
axis([1520 1580 0 4])

m = matlab2tikz(gca);
m.write_tables('sim_se', 'same x')
            
% Sin = S.nlin_sfn.S;
% Sin.P = Sin.P + S.nlin_sfn.num.Pase;
% ASEb = Channels(Signal.wavelength, 0, 'backward'); % zero backward ASE at every iteration
% GaindB = E.propagate(Pump, Sin, ASEb, ASEb, df, 'two-level', 50, false);

save sim_Pp=35mW
            