function WaveForms = get_dso_waveform(N, filepath)

if not(exist('N', 'var'))
    N = 2^19;
end
    
% Settings for DSO
AgilentScope.IPAddr     = '192.168.7.236';      % IP of the Agilent oscilloscope
AgilentScope.WfmLength  = N; % 2^17;                 % Samples per waveform
AgilentScope.SamplingRate = 80E9;               % Samples / second


WaveForms = getAgilentWaveform( AgilentScope ); % Acquire and download


if exist('filepath', 'var')
    disp('Saving waveforms to file:')
    filepath
    save(filepath, 'WaveForms')
end