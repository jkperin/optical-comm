%% Capture BER vs OSNR curve
clear, clc, close all

addpath f/
addpath ../f/
addpath ../mpam/

%%% Experiment parameters
filename = 'PAM4_56G_BER_vs_OSNR_0km_predist_set1';
foldername = 'data/waveforms/BER_vs_OSNR_km/';
dacfile = 'data/waveforms/pam4_rect_Rb=56Gbps_preemph_predist.mat';
S = load(dacfile);
mpam = S.mpam;

% Settings for DSO
AgilentScope.IPAddr     = '192.168.7.236';      % IP of the Agilent oscilloscope
AgilentScope.WfmLength  = 2^21; % 2^17;                 % Samples per waveform
AgilentScope.SamplingRate = 80E9;               % Samples / second

v=visa('agilent', ['TCPIP0::' AgilentScope.IPAddr '::inst0::INSTR']);

captures = 1;
figure(1), box on, hold on
xlabel('OSNR (dB)')
ylabel('log_{10}(BER)')
while true
    fprintf('-- Capture #%d\n', captures)
    
    try
        OSNRdBk = input('OSNR in dB (enter -1 to exit): ');
    catch e
        disp('Invalid OSNR value')
        continue
    end
    
    if OSNRdBk == -1
        break;
    else
        OSNRdB(captures) = OSNRdBk;
    end
    
    WaveForms{captures} = getAgilentWaveform( AgilentScope, v ); % Acquire and download
    
    disp('Waveforms captured. Evaluating BER.')
    
    [ber_count(captures), Vset{captures}] = process_pam_waveforms(WaveForms{captures}, dacfile);
        
    V = Vset{captures};
    
    ber_mzmnl(captures) = pam_ber_from_osnr_mzmnl(mpam.M, OSNRdBk, mpam.Rs/2, V(1:2));
    
    figure(1), hold on
    [OSNRsort, idx] = sort(OSNRdB, 'ascend');
    OSNRplot = linspace(OSNRsort(1), OSNRsort(end)+10, 50);
    ber_ideal = pam_ber_from_osnr(mpam.M, OSNRplot, mpam.Rs/2);
    plot(OSNRplot, log10(ber_ideal), '-b', 'LineWidth', 2)
    plot(OSNRsort, log10(ber_mzmnl(idx)), 'or', 'LineWidth', 2)
    plot(OSNRsort, log10(ber_count(idx)), 'ok', 'LineWidth', 2)
    axis([OSNRsort(1) OSNRsort(end)+10 -8 0])
    legend('Sig-spont limit', 'MZM NL', 'Meaured')
    captures = captures + 1;
end
    

if length(OSNRdB) > 1
    [OSNRdBsort, idx] = sort(OSNRdB, 'ascend');
    experiment.mpam = mpam;
    experiment.OSNRdB = OSNRdBsort;
    experiment.BER.count = ber_count(idx);
    experiment.BER.mzmnl = ber_mzmnl(idx);
    experiment.BER.ideal = ber_ideal(idx);
    experiment.dacfile = dacfile;
    experiment.WaveForms = WaveForms(idx);
    experiment.Vset = Vset{idx};
    save([foldername filename], 'experiment')
    
    if foldername(end) ~= '/' || foldername(end) ~= '\'
        foldername = [foldername '/'];
    end
    
    disp(['File saved at' foldername filename])
    
    figure(1)
    saveas(gca, filename, 'png')
    disp('Figure saved in current folder')
end

    