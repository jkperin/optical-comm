function [simTable, txTable, rxTable] = generate_summary(sim, tx, rx)
%% Simulation parameters
disp('Simulation parameters summary')

rows = {'Bit rate'; 'Number of Variables'; 'Oversampling to emulate continuous time';...
    'Target BER'};
Variables = {'Rb'; 'Nsymb'; 'Mct'; 'BERtarget'};
Values = [sim.Rb/1e9; sim.Nsymb; sim.Mct; sim.BERtarget];
Units = {'Gb/s'; ''; ''; ''};

simTable = table(Variables, Values, Units, 'RowNames', rows)

%% Transmitter paramaters
disp('Transmitter parameters summary')

rows = {'Power swipe (start)'; 'Power swipe (end)'; 'Transmitter wavelength'; 'Laser RIN'; 'Modulator extinction ratio';...
    'Modulator bandwidth'; 'Modulator chirp parameter'};
Variables = {'PtxdBm'; 'PtxdBm';'lamb'; 'RIN'; 'rexdB'; 'modulator.BW'; 'alpha'};
Values = [tx.PtxdBm(1); tx.PtxdBm(end); tx.lamb*1e9; tx.RIN; tx.rexdB; tx.modulator.BW/1e9; tx.alpha];
Units = {'dBm';'dBm';'nm'; 'dB/Hz'; 'dB'; 'GHz'; ''};

txTable = table(Variables, Values, Units, 'RowNames', rows)

%% Receiver paramaters
disp('Receiver parameters summary')

rows = {'Thermal noise one-sided PSD'; 'Optical filter type'; 'Optical filter order';...
    'Optical filter bandwidth';'Electric filter type'; 'Electric filter order'; 'Electric filter bandwidth';...
    'Equalizer type'; 'Equalizer Ntaps'; 'Equalizer adptation rate'};
    
Variables = {'N0'; 'optfilt.type';'optfilt.order'; 'optfilt.fcnorm'; 'elefilt.type';...
    'elefilt.order'; 'elefilt.fcnorm'; 'eq.type'; 'eq.Ntaps'; 'eq.mu'};
Values = {rx.N0; rx.optfilt.type; rx.optfilt.order; rx.optfilt.fcnorm*sim.fs/1e9; rx.elefilt.type;...
    rx.elefilt.order; rx.elefilt.fcnorm*sim.fs/1e9; rx.eq.type; rx.eq.Ntaps; rx.eq.mu};
Units = {'W/Hz'; ''; ''; 'GHz'; ''; ''; 'GHz'; ''; ''; ''};

rxTable = table(Variables, Values, Units, 'RowNames', rows)