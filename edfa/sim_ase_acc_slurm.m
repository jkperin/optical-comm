function sim_ase_acc_slurm(filename, GFF_period, partial_GFF)
%% Simulate ASE accumulation with perfect GFF only after "GFF_period" spans
% Inputs:
% - filename: name of .mat file containig results from optimization
% - GFF_period: period of ideal GFF correction. If GFF_period = 1, then
% ideal GFF is applied at the end of every span.
% - partial_GFF: boolean. If true, fixed GFF equal to the GFF used in the
% last span is used for simulations

addpath results/
addpath f/
addpath data/
addpath ../f/

verbose = false;

disp(filename)
S = load(filename);

% Other input parameters
GFF_period = round(str2double(GFF_period));
if exist('partial_GFF', 'var') % currently not used
    GFF_P1 = load('sim_ase_acc_Ppump=60mW_GFF_period=1_partial_GFF=0.mat');
    partial_GFF = round(str2double(partial_GFF));
    disp('Using last GFF as default GFF filter')
else 
    partial_GFF = 0;
end
fprintf('========= GFF_period = %d =============\n', GFF_period)

spanAttdB = S.spanAttdB;
OffPower = 1e-12; % power set for OFF channels
df = S.problem.df;  
E = S.nlin_sfn.E;
Pump = S.Pump;
SignalIn = S.nlin_sfn.S.sample(S.nlin_sfn.S.wavelength < 1571e-9);
SignalIn.P(SignalIn.P == 0) = OffPower;
ASEf = Channels(SignalIn.wavelength, 0, 'forward');
nSignal = SignalIn; % signal at the output of the nth span.

nGaindB = zeros(SignalIn.N, S.Nspans); % matrix of signal power at the end of each span 
nSignalPdBm = zeros(SignalIn.N, S.Nspans); % matrix of signal power at the end of each span 
num_gffs = ceil(S.Nspans/GFF_period);
nGFF = zeros(SignalIn.N, num_gffs);
nASE = zeros(SignalIn.N, S.Nspans); % matrix of signal power at the end of each span 

gff_cnt = 1;
for n = 1:S.Nspans
    fprintf('n = %d\n', n);
    ASEb = Channels(SignalIn.wavelength, 0, 'backward'); % set backward ASE to zero at every iteration
    
    figure(1), hold on, box on
    plot(nSignal(n).lnm, nSignal(n).PdBm)
    xlabel('Wavelength (nm)')
    ylabel('Input signal power (dBm)')
    axis([1520 1580 -30 -10])
    
    [GaindB, ~, Psignal_out, Pase, sol] = E.propagate(Pump, nSignal(n), ASEf, ASEb, df, 'two-level', 50, false);
    nSignal(n+1) = nSignal(n);
    nSignal(n+1).P = Psignal_out;
        
    if any(sol.y(:) < 0)
        warning('EDF/propagate: solution contains negative power')
    end
    
    if mod(n, GFF_period) == 0 || n == S.Nspans 
        FdB = min(SignalIn.PdBm - Watt2dBm(Psignal_out) + spanAttdB, 0);
        nGFF(:, gff_cnt) = FdB.';
        gff_cnt = gff_cnt + 1;
    elseif partial_GFF > 0
        gff_idx = min(partial_GFF*ceil(n/partial_GFF), 287);
        FdB = GFF_P1.nGFF(:, gff_idx).';
    else
        FdB = zeros(size(GaindB));
    end
    % turn off channels with gain below attenuation
    % nSignal(n+1).P(GaindB < spanAttdB + 0.1) = OffPower;
    fprintf('Number of ON channels = %d / %d\n',  sum(GaindB > spanAttdB), sum(SignalIn.P ~= OffPower))

    nGaindB(:, n) = GaindB.';
    
    nSignal(n+1).PdBm = Watt2dBm(Psignal_out) + FdB - spanAttdB;
    nSignal(n+1).P(nSignal(n+1).P < OffPower) = OffPower;
    ASEf.PdBm = Watt2dBm(Pase) + FdB - spanAttdB;
    
    nSignalPdBm(:, n) = nSignal(n+1).PdBm.';
    nASE(:, n) = ASEf.PdBm.';
     
    if verbose
        figure(2), hold on, box on
        plot(SignalIn.lnm, GaindB - spanAttdB)
        xlabel('Wavelength (nm)')
        ylabel('Gain - Span att. (dB)')
        xlim([1520 1580])

        figure(3), hold on, box on
        plot(SignalIn.lnm, FdB)
        xlabel('Wavelength (nm)')
        ylabel('GFF (dB)')
        xlim([1520 1580])

        figure(4), hold on, box on
        plot(SignalIn.lnm, Watt2dBm(Psignal_out))
        xlabel('Wavelength (nm)')
        ylabel('Output signal power (dBm)')
        axis([1520 1580 -20 0])

        figure(5), hold on, box on
        plot(SignalIn.lnm, ASEf.PdBm)
        xlabel('Wavelength (nm)')
        ylabel('ASE power (dBm)')
        xlim([1520 1580])
        
        drawnow
    end
        
    fprintf('Total input power = %.2f dBm\n', Watt2dBm(sum(nSignal(n+1).P)))
    fprintf('Total ASE power = %.2f dBm\n', Watt2dBm(sum(ASEf.P)))
end

SNR = nSignal(n+1).P./(ASEf.P + S.nlin_sfn.num.NL);
SE = 2*log2(1 + S.problem.Gap*SNR);

fprintf('Total capacity = %.2f Tb/s | %.2f Tb/s\n', sum(S.nlin_sfn.num.SE)*df/1e12, sum(SE)*df/1e12)

if verbose 
    figure, box on, hold on
    plot(SignalIn.lnm, SE)
    plot(SignalIn.lnm, S.nlin_sfn.num.SE)
    plot(SignalIn.lnm, S.nlin_sfn.approx.SE)
    xlabel('Wavelength (nm)')
    ylabel('Spectrale efficiency (bit/s/Hz)')
    legend('Simulation', 'Predicted', 'Predicted approx')

end

output_filename = sprintf('sim_ase_acc_Ppump=%dmW_GFF_period=%d_partial_GFF=%d', round(Pump.P*1e3), GFF_period, partial_GFF);
disp(output_filename)
save(output_filename)
  

% Sin = S.nlin_sfn.S;
% Sin.P = Sin.P + S.nlin_sfn.num.Pase;
% ASEb = Channels(Signal.wavelength, 0, 'backward'); % zero backward ASE at every iteration
% GaindB = E.propagate(Pump, Sin, ASEb, ASEb, df, 'two-level', 50, false);
            
