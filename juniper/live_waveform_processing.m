%% Get DSO data constantly and plot signal after equalization
clear, clc

addpath f/
addpath ../f/
addpath ../mpam/

% Settings for DSO
AgilentScope.IPAddr     = '192.168.7.236';      % IP of the Agilent oscilloscope
AgilentScope.WfmLength  = 2^18; % 2^17;                 % Samples per waveform
AgilentScope.SamplingRate = 80E9;               % Samples / second

v=visa('agilent', ['TCPIP0::' AgilentScope.IPAddr '::inst0::INSTR']);

% dacfile = 'data\waveforms\BER_vs_OSNR_km\pam4_rect_Rb=55Gbps_preemph.mat';
% dacfile = 'data/waveforms/pam4_rect_Rb=56Gbps_preemph_duobin';
dacfile = 'data/waveforms/pam4_rect_Rb=56Gbps_preemph_predist'
waitingTime = 0; % s

Dac = load(dacfile);
OSNRdB = 30;
ber_gauss = -log10(pam_ber_from_osnr(Dac.mpam.M, OSNRdB, Dac.mpam.Rs/2));

iterations = 1;
figure(1), clf, box on, grid on
xlabel('Iterations')
ylabel('-log_{10}(BER)')
while true
    WaveForms = getAgilentWaveform( AgilentScope, v ); % Acquire and download
    
    [ber_count(iterations), Vset] = process_pam_waveforms(WaveForms, dacfile);
           
    figure(1), hold on
    stem(1:iterations, -log10(ber_count), 'b')
    plot([1 iterations+10], ber_gauss*[1 1], ':k')
    axis([1 iterations+10 0 ceil(ber_gauss)])
    legend('Counted', 'Sig-spont limit', 'Location', 'NorthWest')
    iterations = iterations + 1;
    pause(waitingTime)
end
    