clear, clc, close all

addpath f/
addpath ../f/
addpath ../ofdm/
addpath ../ofdm/f/
addpath ../soa/

%% DAC
Tx.DAC.fs = 2687.5*32*1e6; % DAC sampling rate
Tx.DAC.resolution = 8; % DAC effective resolution in bits
Tx.DAC.Vswing = 1; % controls clipping. If Vswing > 1 signal will be clipped
% Tx.DAC.filt = design_filter('butter', 3, 20e9/(sim.fs/2)); % DAC analog frequency response

Rb_tentative = 43e9*1.3;
sim.M = 16;
Rs_tentative = 2*Rb_tentative/log2(sim.M);
Ncp = 14;
Nc = 256;
Nu = 208;
fs_tentative = Rs_tentative*(Nc + Ncp)/Nu;
ofdm_ros = Tx.DAC.fs/fs_tentative;
[p, q] = rat(ofdm_ros)

%% Simulation parameters
sim.Rb = Rb_tentative;
sim.Nsymb = 204; % Number of symbols in pattern
sim.Mct = 1;      
sim.BERtarget = 1.8e-4;
sim.qunatiz = true;
sim.overwrite = true;

N_expected = ceil((Ncp+Nc)*sim.Nsymb*ofdm_ros);
Npad_expected = 2^(ceil(log2(N_expected))) - N_expected
logN = log2(N_expected + Npad_expected)

%% OFDM 
% OFDM constructor: ofdm(Nc, Nu, CS, Rb)
% Nc : number of subcarriers
% Nu : number of used subcarriers
% CS : nominal constellation size
% Rb : bit rate (b/s)
ofdm = ofdm(Nc, Nu, sim.M, sim.Rb, 'palloc'); 
ofdm.set_cyclic_prefix(Ncp/2, Ncp/2);

%% Power allocation 
Hpreemph = 10.^(polyval([-0.0013 0.5846 1.5859], ofdm.fc/1e9)/20); 
Hmod = 1./Hpreemph; % MZM frequency response
Fiber = fiber(10e3);
Himdd = Fiber.Himdd(ofdm.fc, 1550e-9, 0, 'large signal');

PrxdBm = 13; % expected received power in dBm
OSNRdB = 30;
BWref = 12.5e9; % reference bandwidth used to measure OSNR
OSNR = 10^(OSNRdB/10);
N0 = dBm2Watt(PrxdBm)/(BWref*OSNR); % Amplifier one-sided ASE PSD 
varSigSpont = 2*dBm2Watt(PrxdBm)*N0*ofdm.fs/2;

ofdm.power_allocation(Hmod.*Himdd, varSigSpont*ones(size(ofdm.fc))/ofdm.Nc, sim.BERtarget, true);

% Try detecting
AdEq.mu = 1e-3;
AdEq.Ntrain = Inf; % Number of symbols used in training (if Inf all symbols are used)
[xd, AdEq.trainSeq] = ofdm.signal(sim.Nsymb); 
Xn = ofdm.detect(xd, AdEq, true);
figure, plot(ofdm.fc/1e9, abs(Xn).^2)
ofdm.countBER([1 1], false) % sanity check

%% Resample to DAC sampling rate
xkDACraw = resample(xd, p, q);

%% Zero pad
Nxd = length(xkDACraw);
Npad = 2^(ceil(log2(Nxd)))-Nxd;
fprintf('Generated sequence has %d points. Padding with %d zeros to get to next power of 2\n', Nxd, Npad)
xkDAC = [zeros(1, floor(Npad/2)) xkDACraw zeros(1, ceil(Npad/2))];

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

figure, plot(xq)

%% Save file
filename = sprintf('data/waveforms/ofdm%d_Rb=%dGbps_%d', ofdm.CS, round(ofdm.Rb/1e9), sim.Nsymb);

filename
if exist([filename '.mat'], 'file') && not(sim.overwrite)
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
    fsrfDataJess = JESS_RAMFill_General(fsrfDataJess, xq, zeros(size(xq)), zeros(size(xq)), zeros(size(xq)));
end
