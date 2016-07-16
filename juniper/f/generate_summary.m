function generate_summary(mpam, Tx, Fibers, EDFA, Rx, sim)
%% Simulation parameters
disp('Simulation parameters summary')

rows = {'Bit rate'; 'Number of symbols'; 'Oversampling to emulate continuous time';...
    'Oversampling ratio of transmitter DSP'; 'Oversampling ratio of receiver DSP'; 'Target BER'};
Variables = {'Rb'; 'Nsymb'; 'Mct'; 'ros.txDSP'; 'ros.rxDSP'; 'BERtarget'};
Values = [sim.Rb/1e9; sim.Nsymb; sim.Mct; sim.ros.txDSP; sim.ros.rxDSP; sim.BERtarget];
Units = {'Gb/s'; ''; ''; ''; ''; ''};

simTable = table(Variables, Values, Units, 'RowNames', rows)

%% M-PAM
mpam.summary();

%% Transmitter paramaters
disp('Transmitter parameters summary')

rows = {'TX DSP oversampling ratio'; 'DAC sampling ratio'; 'DAC effective resolution';... 
    'Modulator type'; 'Modulator bandwidth'; 'Modulator extinction ratio'; 'Modulator chirp parameter'};
Variables = {'DAC.ros'; 'DAC.fs';'DAC.resolution'; 'Mod.type'; 'Mod.BW'; 'rexdB'; 'alpha'};
Values = {Tx.DAC.ros; Tx.DAC.fs/1e9; Tx.DAC.resolution; Tx.Mod.type; Tx.Mod.BW/1e9; Tx.rexdB; Tx.alpha};
Units = {''; 'GS/s';'bits'; ''; 'GHz'; 'dBm'; ''};

txTable = table(Variables, Values, Units, 'RowNames', rows)

%% Laser
Tx.Laser.summary();

%% Fibers
for k = 1:length(Fibers)
    Fibers(k).summary(Tx.Laser.lambda);
end

%% EDFA
for k = 1:length(EDFA)
    EDFA(1).summary();
end

%% Receiver paramaters
disp('Receiver parameters summary')

rows = {'RX DSP oversampling ratio'; 'Thermal noise one-sided PSD'; 'Optical filter type';...
    'Optical filter bandwidth';...
    'Equalizer type'; 'Equalizer Ntaps'; 'Equalizer adptation rate'};
    
Variables = {'ADC.ros'; 'N0'; 'optfilt.type'; 'optfilt.fcnorm'; 'eq.type'; 'eq.Ntaps'; 'eq.mu'};
Values = {Rx.ADC.ros; Rx.N0; Rx.optfilt.type; Rx.optfilt.fcnorm*sim.fs/1e9; Rx.eq.type; Rx.eq.Ntaps; Rx.eq.mu};
Units = {''; 'W/Hz'; ''; 'GHz'; ''; ''; ''};

rxTable = table(Variables, Values, Units, 'RowNames', rows)

%% Photodiode
Rx.PD.summary();

%% ADC
if isfield(Rx, 'ADC')
    disp('-- ADC parameters summary:')

    rows = {'Quantization on?'; 'Effective resolution'; 'Sampling rate'; 'Antialiasing filter type';...
        'Antialiasing filter order'; 'Antialiasing filter bandwidth'; 'Clipping ratio'};

    Variables = {'sim.quantiz'; 'ENOB'; 'fs'; 'filt.type';...
        'filt.order'; 'filt.fcnorm'; 'rclip'};
    
    Values = {sim.quantiz; Rx.ADC.ENOB; Rx.ADC.fs/1e9; Rx.ADC.filt.type;...
        Rx.ADC.filt.order; Rx.ADC.filt.fcnorm*sim.fs/1e9; Rx.ADC.rclip};
    Units = {''; 'bits'; 'GS/s'; ''; ''; 'GHz'; '%'};

    ADCTable = table(Variables, Values, Units, 'RowNames', rows)
end