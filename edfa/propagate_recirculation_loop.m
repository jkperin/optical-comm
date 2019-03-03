%% Simulate propagation in recirculation loop
clear, clc, close all

addpath results/
addpath f/
addpath data/
addpath ../f/

verbose = true;

filename = 'results/capacity_vs_pump_power_EDF=corning high NA_pump=60mW_980nm_ChDf=50GHz_L=264_x_50km.mat';
S = load(filename);
GFF_P1 = load('sim_ase_acc_Ppump=60mW_GFF_period=1_partial_GFF=0.mat');

N_spans_per_loop = 5;
N_loops = floor(220 / N_spans_per_loop);
AdditionalLossPerLoopdB = 6.5 + 3.5 + 5; % AOM & 3dB coupler + Polarization scrambler + Waveshaper
AdditionalAmpGaindB = AdditionalLossPerLoopdB; % gain of additional amp compensates for additional losses
AdditionalAmpNF = 5; % noise figure of additional amp in loop
% Noise of additional amplifier
% 2*nsp*(G-1)/G*hv*df: Note that this noise is already attenuated by 1/G
AdditionalNoise = 10^(AdditionalAmpNF/10)*S.Signal.Ephoton*S.problem.df; 

spanAttdB = S.spanAttdB;
OffPower = 1e-12; % power set for OFF channels

E = S.lin.E;
Pump = S.Pump;
SignalIn = S.nlin.S;
SignalIn.P(SignalIn.P == 0) = OffPower;
OptPowerProfiledBm = SignalIn.PdBm; % optimal power profile in dBm
Signal = SignalIn;

ASEf = Channels(SignalIn.wavelength, AdditionalNoise, 'forward');
ASEf.df = S.problem.df;


FixGFFperAmp = GFF_P1.nGFF(:, end).'; % make fix GFF equal to last GFF in loop

nGaindB = zeros(SignalIn.N, S.Nspans); % matrix of signal power at the end of each span 
nSignalPdBm = zeros(SignalIn.N, S.Nspans); % matrix of signal power at the end of each span 
nASE = zeros(SignalIn.N, S.Nspans); % matrix of signal power at the end of each span 
nWSS = zeros(SignalIn.N, N_loops); 
lamb_lim = 1e9 * SignalIn.wavelength([1 end]);

span = 1;
for loop = 1:N_loops
    for k = 1:N_spans_per_loop
        fprintf('loop %d x span %d = %d spans\n', loop, k, span);
         
        Signal.plot(figure(1));
        axis([lamb_lim -30 -10])
        ylabel('Input signal power (dBm)')
        title(sprintf('span = %d', span))
        
        %%
        [Signal, ASEf, GaindB] = Amp_GFF_Loss(E, Pump, Signal, ASEf, spanAttdB,...
            OptPowerProfiledBm, FixGFFperAmp);
     
        % log gain, power, and ASE after each span
        nGaindB(:, span) = GaindB.';
        nSignalPdBm(:, span) = Signal.PdBm.';
        nASE(:, span) = ASEf.PdBm.';
        
        % Plot
        if verbose
            figure(2), hold on, box on
            plot(SignalIn.lnm, GaindB - spanAttdB)
            xlabel('Wavelength (nm)')
            ylabel('Gain - Span att. (dB)')
            xlim(lamb_lim)
            title(sprintf('span = %d', span))

            figure(3), hold on, box on
            plot(SignalIn.lnm, ASEf.PdBm)
            xlabel('Wavelength (nm)')
            ylabel('ASE power (dBm)')
            xlim(lamb_lim)
            title(sprintf('span = %d', span))
            
            drawnow
        end
        
        span = span + 1;
    end
    
    %% Additional amplifier and span
%     [Signal, ASEf, GaindB, WSS] = Amp_GFF_Loss(E, Pump, Signal, ASEf, spanAttdB,...
%         OptPowerProfiledBm);
% 
%     % log gain, power, and ASE after each span
%     nGaindB(:, span) = GaindB.';
%     nSignalPdBm(:, span) = Signal.PdBm.';
%     nASE(:, span) = ASEf.PdBm.';
%     nWSS(:, loop) = WSS.';
%     
%     span = span + 1;
    
    %% Amp and loss
    Signal.PdBm = Signal.PdBm - AdditionalLossPerLoopdB + AdditionalAmpGaindB;
    
    % Account for noise figure of additional amplifier    
    ASEf.PdBm = ASEf.PdBm - AdditionalLossPerLoopdB + AdditionalAmpGaindB;
    ASEf.P = ASEf.P + AdditionalNoise;
    
    FdB = min(OptPowerProfiledBm - Signal.PdBm, 0);   
    Signal.PdBm = Signal.PdBm + FdB;
    ASEf.PdBm = ASEf.PdBm + FdB;
    nWSS(:, loop) = FdB.';
    
    %% 
    if verbose
        figure(4), hold on, box on
        plot(SignalIn.lnm, nWSS(:, loop))
        xlabel('Wavelength (nm)')
        ylabel('WSS gain profile (dB)')
        xlim(lamb_lim)
        title(sprintf('span = %d', span))
    end
    
    1;
end

%%
SNR = Signal.P./(ASEf.P + S.nlin.num.NL);
SE = 2*log2(1 + S.problem.Gap*SNR);

df = S.problem.df;
fprintf('Total capacity = %.2f Tb/s | %.2f Tb/s\n', sum(S.nlin.num.SE)*df/1e12, sum(SE)*df/1e12)

if verbose 
    figure, box on, hold on
    plot(SignalIn.lnm, SE)
    plot(SignalIn.lnm, S.nlin.num.SE)
    plot(SignalIn.lnm, S.nlin.approx.SE)
    xlabel('Wavelength (nm)')
    ylabel('Spectrale efficiency (bit/s/Hz)')
    legend('Simulation', 'Predicted', 'Predicted approx')
end

% output_filename = sprintf('sim_ase_acc_Ppump=%dmW_GFF_period=%d_partial_GFF=%d', round(Pump.P*1e3), GFF_period, partial_GFF);
% disp(output_filename)
% save(output_filename)