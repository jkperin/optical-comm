function ofdm_simulation_summary(sim, ofdm, tx, fibers, rx)
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

rows = {'Bit rate'; 'OFDM type'; 'Number of symbols'; 'Oversampling to emulate continuous time';...
    'Target BER'; 'Include phase noise?'; 'Inlcude RIN?'; 'Include PMD?'; 'Include quantization?'};

Variables = {'Rb'; 'OFDM'; 'Nsymb'; 'Mct';...
    'BERtarget'; 'phase_noise'; 'RIN'; 'PMD'; 'quantiz'};

Values = {sim.Rb/1e9; sim.OFDM; sim.Nsymb; sim.Mct;...
    sim.BERtarget; sim.phase_noise; sim.RIN; sim.PMD; sim.quantiz};

Units = {'Gb/s'; '';...
    ''; '';...
    ''; ''; ''; ''; ''};

simTable = table(Variables, Values, Units, 'RowNames', rows)

%% Modulation format
ofdm.summary()

%% Transmitter laser parameters
tx.Laser.summary()

%% Transmitter paramaters
disp('-- Transmitter parameters summary:')

rows = {'Modulator type'; 'Modulator bandwidth'; 'Transmitter filter type';...
    'Transmitter filter order'; 'Transmitter filter bandwidth'};

Variables = {'Mod.type'; 'Mod.BW'; 'filt.type'; 'filt.order'; 'filt.fcnorm'};

Values = {tx.Mod.type; tx.Mod.BW/1e9; tx.Mod.filt.type; tx.Mod.filt.order; tx.Mod.filt.fcnorm*sim.fs/2e9};
Units = {''; 'GHz'; ''; ''; 'GHz'};

txTable = table(Variables, Values, Units, 'RowNames', rows)

%% Fiber
for k = 1:length(fibers)
    fprintf('-- Fiber #%d\n', k)
    fibers(k).summary(tx.Laser.lambda)
end

%% Receiver paramaters
disp('-- Receiver parameters summary')

rows = {'Thermal noise one-sided PSD'};
Variables = {'N0'};
Values = [rx.N0];
Units = {'W/Hz';};

rxTable = table(Variables, Values, Units, 'RowNames', rows)

%% Photodiode
rx.PD.summary()

%% ADC
if isfield(rx, 'ADC')
    disp('-- ADC parameters summary:')

    rows = {'Quantization on?'; 'Effective resolution'; 'Sampling rate'; 'Antialiasing filter type';...
        'Antialiasing filter order'; 'Antialiasing filter bandwidth'};

    Variables = {'sim.quantiz'; 'ENOB'; 'fs'; 'filt.type';...
        'filt.order'; 'filt.fcnorm'};
    
    Values = {sim.quantiz; rx.ADC.ENOB; rx.ADC.fs/1e9; rx.ADC.filt.type;...
        rx.ADC.filt.order; rx.ADC.filt.fcnorm*sim.fs/2e9};
    Units = {''; 'bits'; 'GS/s'; ''; ''; 'GHz'};

    ADCTable = table(Variables, Values, Units, 'RowNames', rows)
end

%% Adaptive Equalization
if isfield(rx, 'AdEq')
    disp('-- Adaptive equalization summary:')

    rows = {'Training sequence length'; 'Adaptation rate'};

    Variables = {'Ntrain'; 'mu'};
    
    Values = {rx.AdEq.Ntrain; rx.AdEq.mu};
    
    Units = {''; ''};

    AdEqTable = table(Variables, Values, Units, 'RowNames', rows)
end