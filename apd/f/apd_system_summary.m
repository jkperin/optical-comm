function apd_system_summary(mpam, Tx, Fiber, Rx, Apd, sim)

%% Simulation parameters
disp('Simulation parameters summary')

rows = {'Bit rate'; 'Number of symbols'; 'Oversampling to emulate continuous time';...
    'Oversampling ratio of transmitter DSP'; 'Oversampling ratio of receiver DSP'; 'Target BER'};
Variables = {'Rb'; 'Nsymb'; 'Mct'; 'ros.txDSP'; 'ros.rxDSP'; 'BERtarget'};
Values = [sim.Rb/1e9; sim.Nsymb; sim.Mct; sim.ros.txDSP; sim.ros.rxDSP; sim.BERtarget];
Units = {'Gb/s'; ''; ''; ''; ''; ''};

simTable = table(Variables, Values, Units, 'RowNames', rows)

%% M-PAM
mpam.summary()

%% Transmitter paramaters
disp('Transmitter parameters summary')

rows = {'TX DSP oversampling ratio'; 'DAC sampling ratio'; 'DAC effective resolution';... 
    'Modulator type'; 'Modulator bandwidth'; 'Modulator extinction ratio'; 'Modulator chirp parameter'};
Variables = {'DAC.ros'; 'DAC.fs';'DAC.resolution'; 'Mod.type'; 'Mod.BW'; 'rexdB'; 'alpha'};
Values = {Tx.DAC.ros; Tx.DAC.fs/1e9; Tx.DAC.resolution; Tx.Mod.type; Tx.Mod.BW/1e9; Tx.rexdB; Tx.alpha};
Units = {''; 'GS/s';'bits'; ''; 'GHz'; 'dBm'; ''};

txTable = table(Variables, Values, Units, 'RowNames', rows)

%% Transmitter laser
Tx.Laser.summary()

%% Fiber
Fiber.summary()

%% Receiver
%% Receiver paramaters
disp('Receiver parameters summary')

rows = {'RX DSP oversampling ratio'; 'Thermal noise one-sided PSD'; 'Receiver filtering';...
    'Equalizer type'; 'Equalizer Ntaps'};
    
Variables = {'ADC.ros'; 'N0'; 'filtering'; 'eq.type'; 'eq.Ntaps'};
Values = {Rx.ADC.ros; Rx.N0; Rx.filtering; Rx.eq.type; Rx.eq.Ntaps};
Units = {''; 'W/Hz'; ''; ''; ''};

rxTable = table(Variables, Values, Units, 'RowNames', rows)

%% APD
Apd.summary()