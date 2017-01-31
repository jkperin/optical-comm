function coherent_simulation_summary(sim, tx, fiber, rx)
warning('off', 'MATLAB:table:DuplicateDimnamesVarnamesBackCompat');

%% Simulation parameters
disp('-- Simulation parameters summary:')

if ~isfield(sim, 'RIN')
    sim.RIN = false;
end

if ~isfield(sim, 'quantiz')
    sim.quantiz = false;
end

if ~isfield(sim, 'PMD')
    sim.PMD = false;
end

if ~isfield(sim, 'phase_noise')
    sim.phase_noise = false;
end

rows = {'Bit rate'; 'Number of polarizations'; 'Modulation format'; 'Modulation Order';...
    'Symbol Rate'; 'Number of symbols'; 'Oversampling to emulate continuous time';...
    'Target BER'; 'Include phase noise?'; 'Inlcude RIN?'; 'Include PMD?'; 'Include quantization?'};

Variables = {'Rb'; 'Npol'; 'ModFormat.type'; 'ModFormat.M';...
    'ModFormat.Rs'; 'Nsymb'; 'Mct';...
    'BERtarget'; 'phase_noise'; 'RIN'; 'PMD'; 'quantiz'};

Values = {sim.Rb/1e9; sim.Npol; sim.ModFormat.type; sim.ModFormat.M;  
    sim.ModFormat.Rs/1e9; sim.Nsymb; sim.Mct;...
    sim.BERtarget; sim.phase_noise; sim.RIN; sim.PMD; sim.quantiz};

Units = {'Gb/s'; ''; ''; '';...
    'Gbaud'; ''; '';...
    ''; ''; ''; ''; ''};

simTable = table(Variables, Values, Units, 'RowNames', rows)

%% Modulation format
sim.ModFormat.summary()

%% Transmitter laser parameters
tx.Laser.summary()

%% Transmitter paramaters
disp('-- Transmitter parameters summary:')

rows = {'Modulator type'; 'Modulator bandwidth'; 'Transmitter filter type';...
    'Transmitter filter order'; 'Transmitter filter bandwidth'};

Variables = {'Mod'; 'Mod.BW'; 'filt.type'; 'filt.order'; 'filt.fcnorm'};

Values = {sim.Modulator; tx.Mod.BW/1e9; tx.filt.type; tx.filt.order; tx.filt.fcnorm*sim.fs/1e9};
Units = {''; 'GHz'; ''; ''; 'GHz'};

txTable = table(Variables, Values, Units, 'RowNames', rows)

%% Fiber
fiber.summary(tx.Laser.lambda)

%% Receiver paramaters
disp('-- Receiver parameters summary')

rows = {'Thermal noise one-sided PSD'};
Variables = {'N0'};
Values = [rx.N0];
Units = {'W/Hz';};

rxTable = table(Variables, Values, Units, 'RowNames', rows)

%% Optical amplifier
if isfield(rx, 'OptAmp')
    disp('-- Simulation includes optical amplifier');
    rx.OptAmp.summary()
end

%% Local Oscillator
rx.LO.summary()

%% Photodiode
rx.PD.summary()

%% ADC
if isfield(rx, 'ADC')
    disp('-- ADC parameters summary:')

    rows = {'Quantization on?'; 'Effective resolution'; 'Sampling rate'; 'Antialiasing filter type';...
        'Antialiasing filter order'; 'Antialiasing filter bandwidth'; 'Clipping ratio'};

    Variables = {'sim.quantiz'; 'ENOB'; 'fs'; 'filt.type';...
        'filt.order'; 'filt.fcnorm'; 'rclip'};
    
    Values = {sim.quantiz; rx.ADC.ENOB; rx.ADC.fs/1e9; rx.ADC.filt.type;...
        rx.ADC.filt.order; rx.ADC.filt.fcnorm*sim.fs/1e9; rx.ADC.rclip};
    Units = {''; 'bits'; 'GS/s'; ''; ''; 'GHz'; '%'};

    ADCTable = table(Variables, Values, Units, 'RowNames', rows)
end

%% Adaptive Equalization
if isfield(rx, 'AdEq')
    disp('-- Adaptive equalization summary:')

    rows = {'Type'; 'Structure'; 'Number of taps'; 'Oversampling ratio'; 
        'Training sequence length'; 'Adaptation rate'};

    Variables = {'type'; 'structure'; 'Ntaps'; 'ros';...
        'Ntrain'; 'mu'};
    
    Values = {rx.AdEq.type; rx.AdEq.structure; rx.AdEq.Ntaps; rx.AdEq.ros;
        rx.AdEq.Ntrain; rx.AdEq.mu};
    
    Units = {''; ''; ''; ''; ''; ''};

    AdEqTable = table(Variables, Values, Units, 'RowNames', rows)
end

%% Carrier phase recovery
if isfield(rx, 'CPR')
    disp('-- Carrier phase recovery summary:')

    if strcmpi(rx.CPR.type, 'feedforward')
        if strcmpi(rx.CPR.phaseEstimation, 'NDA')
            phaseEstimation = [rx.CPR.phaseEstimation '-' rx.CPR.NDAorder];
        else
            phaseEstimation = rx.CPR.phaseEstimation;
        end
        rows = {'Type'; 'Phase Estimation'; 'Delay'; 'Set up time'; 
            'Filter'; 'FIR/IIR'; 'Structure'; 'Number of taps'};

        Variables = {'type'; 'phaseEstimation'; 'Delay'; 'Ntrain'; 'Filter';...
            'FilterType'; 'structure'; 'Ntaps'};
        
        Values = {rx.CPR.type; phaseEstimation; rx.CPR.Delay; rx.CPR.Ntrain;...
            rx.CPR.Filter; rx.CPR.FilterType; rx.CPR.structure; rx.CPR.Ntaps};
        Units = {''; ''; 'Symbols'; 'Symbols'; ''; ''; ''; ''};
    else
        rows = {'Type'; 'Phase Estimation'; 'Delay'; 'Set up time'; 
            'CT to DT conversion method'; 'Loop filter damping'; 'Loop filter relaxataion frequency'};

        Variables = {'type'; 'phaseEstimation'; 'Delay'; 'Ntrain';...
            'CT2DT'; 'csi'; 'wn'};
        Values = {rx.CPR.type; rx.CPR.phaseEstimation; rx.CPR.Delay; rx.CPR.Ntrain;...
            rx.CPR.CT2DT; rx.CPR.csi; rx.CPR.wn/1e9};
        Units = {''; ''; 'Symbols'; 'Symbols'; ''; ''; 'Grad/s'};        
    end
        
    CPRTable = table(Variables, Values, Units, 'RowNames', rows)
end

%% Analog
if isfield(rx, 'Analog') && strcmpi(sim.ModFormat.type, 'QAM')
    disp('-- Analog receiver parameters summary')

    rows = {'Antialiasing filter type'; 'Antialiasing filter order'; 'Antialiasing filter bandwidth';...
        'Carrier phase recovery method'; 'Phase estimation method'; 'Delay'; 'Loop filter damping factor'};

    Variables = {'filt.type'; 'filt.order'; 'filt.fcnorm';...
        'CarrierPhaseRecovery'; 'CPRmethod'; 'Additional loop delay'; 'csi'};
    Values = {rx.Analog.filt.type; rx.Analog.filt.order; rx.Analog.filt.fcnorm*sim.fs/1e9;...
        rx.Analog.CarrierPhaseRecovery; rx.Analog.CPRmethod; rx.Analog.Delay*1e12; rx.Analog.csi};

    Units = {''; ''; 'GHz'; ''; ''; 'ps'; ''};

    AnalogTable = table(Variables, Values, Units, 'RowNames', rows)    
end
