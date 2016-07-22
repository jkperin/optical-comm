function load_waveform(file)
addpath f/
addpath ../f/
addpath ../mpam/

S = load(file);

rows = {'DAC sampling rate'; 'Number of samples'; 'Trigger frequency';...
    'Symbol rate'; 'PAM order'; 'Bit rate';
    'Oversampling ratio of transmitter DSP'};
Variables = {'Tx.DAC.fs'; 'Npoints'; 'ftrigger';...
    'mpam.Rs'; 'mpam.M'; 'mpam.Rb';...
    'sim.ros.txDSP'};
Values = [S.Tx.DAC.fs/1e9; S.Npoints; S.ftrigger/1e9;...
    S.mpam.Rs/1e9; S.mpam.M; S.mpam.Rb/1e9;...
    S.sim.ros.txDSP];
Units = {'GS/s'; ''; 'GHz'; 'GBaud'; ''; 'Gbit/s'; ''};

table(Variables, Values, Units, 'RowNames', rows)


%% Load waveform in DAC
global fsrfDataJess

xq = S.xq;
xqtrig = S.xqtrig;

if isempty(fsrfDataJess)
    disp('Waveform not loaded to DAC')
else
    disp('Waveform loaded to DAC')
    fsrfDataJess = JESS_RAMFill_General(fsrfDataJess, xq, xqtrig, zeros(size(xq)), zeros(size(xq)));
end
