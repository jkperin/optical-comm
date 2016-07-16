%% Generate DAC waveform
clear, clc, close all

addpath f/
addpath ../f/
addpath ../mpam/

%% Pre calculations to determine bitrate, trigger frequency, etc
Rs_tentative = 28e9;
Npoints = 2^12; % desired number of points to load to DAC
Tx.DAC.fs = 2687.5*32*1e6; % DAC sampling rate

k2 = 120;
ftrigger = k2*Tx.DAC.fs/Npoints;
k1 = round(Rs_tentative*Npoints/(k2*Tx.DAC.fs));
Rs = k1*ftrigger;

%% Simulation parameters
sim.M = 4; % PAM order
sim.Rb = Rs*log2(sim.M); % bit rate
sim.ros.txDSP = 2; % oversampling ratio of transmitter DSP. Must be integer
sim.Nzero = 5; % set first and last Nzero symbols to 0
sim.qunatiz = true;
sim.preemph = true;
sim.preemphRange = 20e9;
sim.duobinary = ~true;

%% DAC
Tx.DAC.fs = 2687.5*32*1e6; % DAC sampling rate
Tx.DAC.ros = sim.ros.txDSP; % oversampling ratio of transmitter DSP
Tx.DAC.Vswing = 0.9;
Tx.DAC.resolution = 8; % bits

%% M-PAM
pulse_shape = select_pulse_shape('rect', sim.ros.txDSP);
% pulse_shape = select_pulse_shape('rrc', sim.ros.txDSP, 0.25, 4);
mpam = PAM(sim.M, sim.Rb, 'equally-spaced', pulse_shape);

% estimate number of symbols necessary to minimize padding
sim.Nsymb = ceil(Npoints*mpam.Rs/Tx.DAC.fs);

%% Summary
rows = {'DAC sampling rate'; 'Number of samples'; 'Trigger frequency';...
    'Symbol rate'; 'PAM order'; 'Bit rate';
    'Oversampling ratio of transmitter DSP'};
Variables = {'Tx.DAC.fs'; 'Npoints'; 'ftrigger';...
    'mpam.Rs'; 'mpam.M'; 'mpam.Rb';...
    'sim.ros.txDSP'};
Values = [Tx.DAC.fs/1e9; Npoints; ftrigger/1e9;...
    mpam.Rs/1e9; mpam.M; mpam.Rb/1e9;...
    sim.ros.txDSP];
Units = {'GS/s'; ''; 'GHz'; 'GBaud'; ''; 'Gbit/s'; ''};

table(Variables, Values, Units, 'RowNames', rows)

%% Generate data
dataTX = randi([0 mpam.M-1], [1 sim.Nsymb]);
dataTX(1:sim.Nzero) = 0;
dataTX(end-sim.Nzero+1:end) = 0;

%% Duobinary enconding
if isfield(sim, 'duobinary') && sim.duobinary
    mpamdb = mpam.set_levels(0:mpam.M-1, 0.5 + (0:mpam.M-2));    
    
    xd = mpamdb.mod(dataTX); % Modulated PAM signal
    
    xd = duobinary_encoding(xd);
    
    xd_enc = xd;
    
    ximp = upsample(xd_enc, mpam.pulse_shape.sps);
    xk = filter(ones(1, mpam.pulse_shape.sps)/mpam.pulse_shape.sps, 1, ximp);
    xk = xk - mean(xk);
else
    xk = mpam.signal(dataTX); % Generate signal at sim.ros.txDSP rate
    xk = xk - mean(xk); 
end  

%% Preemphasis
if sim.preemph
    f = freq_time(sim.Nsymb*sim.ros.txDSP, mpam.Rs*sim.ros.txDSP);
    femph = abs(f);
    femph(femph >= sim.preemphRange) = 0;
    emphasis_filter = 10.^(polyval([-0.0013 0.5846 1.5859], femph/1e9)/20);    

    figure, plot(f/1e9, 20*log10(abs(emphasis_filter)))
    xlabel('Frequency (GHz)')
    ylabel('Preemphasis filter (dB)')

    xk = real(ifft(fft(xk).*ifftshift(emphasis_filter)));

    figure, eyediagram(xk, 32)
    title('Signal after preemphasis')
end

%% Resample to DAC sampling rate
[p, q] = rat(Tx.DAC.fs/(mpam.Rs*sim.ros.txDSP));
xkDAC = resample(xk, p, q);

%% Generate trigger signal
[~, t] = freq_time(length(xkDAC), Tx.DAC.fs);
xtrig = sin(2*pi*ftrigger*t);
xmin = min(xtrig);
xmax = max(xtrig);
xamp = xmax - xmin;    
dx = xamp/(2^(Tx.DAC.resolution)-1);  
sig_range = xmin:dx:xmax;
partition = sig_range(1:end-1) + dx/2;
[~, xqtrig] = quantiz(0.9*xtrig, partition, 0:255); 

if not(isInteger(log2(length(xkDAC))))
    Np2 = ceil(log2(length(xkDAC)));
    xkDAC = [xkDAC xkDAC(1:Np2-length(xkDAC))];
    xqtrig  = [xktrig xktrig(1:Np2-length(xkDAC))];
    warning('Cyclic extension of %d points was necessary in order to make DAC waveform power of 2', Np2-length(xkDAC))
end

%% Quantization 
xkDAC = xkDAC/sqrt(mean(abs(xkDAC).^2));
xkDAC = xkDAC - mean(xkDAC);
xmin = min(xkDAC);
xmax = max(xkDAC);
xamp = xmax - xmin;    
dx = xamp/(2^(Tx.DAC.resolution)-1);  
sig_range = xmin:dx:xmax;
partition = sig_range(1:end-1) + dx/2;
[~, xq, varQ] = quantiz(Tx.DAC.Vswing*xkDAC, partition, 0:255); 

%% DAC
sim.Mct = 16; % oversampling ratio to emulate continuous time
sim.N = sim.Nsymb*sim.Mct; % number of points 

%% Time and frequency
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'
[sim.f, sim.t] = freq_time(sim.N, sim.fs);

% DAC filter
Tx.DAC.filt = design_filter('bessel', 5, 10e9/(sim.fs/2)); % DAC analog frequency response
xt = dac(xk, Tx.DAC, sim, true);

%% Save file
filename = sprintf('data/waveforms/pam%d_%s_Rb=%dGbps', mpam.M, mpam.pulse_shape.type, round(mpam.Rb/1e9));
if sim.preemph
    filename = [filename '_preemph'];
end

if sim.duobinary
    filename = [filename '_duobin'];
end

filename
if exist([filename '.mat'], 'file')
    disp('File already exists! Using existing file instead')
    load([filename '.mat'])
else
    save(filename)
end

%% Load waveform in DAC
global fsrfDataJess

if isempty(fsrfDataJess)
    disp('Waveform not loaded to DAC')
else
    disp('Waveform loaded to DAC')
    fsrfDataJess = JESS_RAMFill_General(fsrfDataJess, xq, xqtrig, zeros(size(xq)), zeros(size(xq)));
end




