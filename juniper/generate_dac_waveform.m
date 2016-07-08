%% Generate DAC waveform
clear, clc, close all

addpath f/
addpath ../f/
addpath ../mpam/

sim.Nsymb = 2^10;
sim.Mct = 12;
sim.ros.txDSP = 4;
sim.N = sim.Nsymb*sim.Mct;
sim.quantiz = false; % DAC quantization
sim.Ndiscard = 128;

%% DAC
Tx.DAC.fs = 2687.5*32*1e6; % DAC sampling rate
Tx.DAC.ros = sim.ros.txDSP; % oversampling ratio of transmitter DSP
Tx.DAC.resolution = Inf; % DAC effective resolution in bits
% Tx.DAC.filt = design_filter('butter', 3, 20e9/(sim.fs/2)); % DAC analog frequency response
%% 
sim.M = 4;
sim.Rs = Tx.DAC.fs/Tx.DAC.ros;
sim.Rb = log2(sim.M)*sim.Rs;

%% Time and frequency
sim.fs = sim.Rs*sim.Mct;  % sampling frequency in 'continuous-time'

dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt);
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df);

sim.t = t;
sim.f = f;

Tx.DAC.filt = design_filter('butter', 5, 15e9/(sim.fs/2)); % DAC analog frequency response

%% Pulse shape
% pulse_shape = select_pulse_shape('rect', sim.ros.txDSP);
pulse_shape = select_pulse_shape('rrc', sim.ros.txDSP, 0.25, 4);

%% M-PAM
mpam = PAM(sim.M, sim.Rb, 'equally-spaced', pulse_shape);

dataTX = randi([0 mpam.M-1], [1 sim.Nsymb]);

xk = mpam.signal(dataTX);
xk = xk - mean(xk);

%% Duobinary enconding
if isfield(sim, 'duobinary') && sim.duobinary
    mpamdb = mpam.set_levels(0:mpam.M-1, 0.5 + (0:mpam.M-2));    
    
    xd = mpamdb.signal(dataTX); % Modulated PAM signal
    
    xd = duobinary_encoding(xd);
    
    xd_enc = xd;
else
    % Ajust levels to desired transmitted power and extinction ratio
    mpam = mpam.adjust_levels(Tx.Ptx, Tx.rexdB);
    mpam = mpam.norm_levels(); % presevers extinction ratio
    
    xd = mpam.signal(dataTX); % Modulated PAM signal
end  

%% Preemphasis


%% Quantization 
xmin = min(xk);
xmax = max(xk);
Vswing = 0.8;

xq = round(256*Vswing*xk/(xmax-xmin));

save_waveform(xq, sprintf('pam4_%s_waveform.txt', mpam.pulse_shape.type))

%% DAC
xt = dac(xk, Tx.DAC, sim, true);




