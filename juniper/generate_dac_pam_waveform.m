%% Generate DAC waveform
clear, clc, close all

addpath f/
addpath ../f/
addpath ../mpam/

% S = load('data/waveforms/BER_vs_OSNR_b2b_MZM_predist/pam4_rect_Rb=55Gbps_preemph_predist.mat');

%% Pre calculations to determine bitrate, trigger frequency, etc
Rs_tentative = 28e9;
Npoints = 3*2^12; % desired number of points to load to DAC
Tx.DAC.fs = 2625*32*1e6; % DAC sampling rate 2687.5

k2 = 128; % number of trigger cycles within DAC waveform window
ftrigger = k2*Tx.DAC.fs/Npoints;
k1 = round(Rs_tentative/ftrigger);
Rs = k1*ftrigger

%% Simulation parameters
sim.M = 4; % PAM order
sim.Rb = Rs*log2(sim.M); % bit rate
sim.ros.txDSP = 3; % oversampling ratio of transmitter DSP. Must be integer
sim.Nzero = 5; % set first and last Nzero symbols to 0
sim.qunatiz = true;
sim.preemph = true;
sim.preemphRange = 25e9;
sim.mzm_predistortion = true;
sim.duobinary = false;
sim.overwrite = true;

%% DAC
Tx.DAC.ros = sim.ros.txDSP; % oversampling ratio of transmitter DSP
Tx.DAC.Vswing = 1;
Tx.DAC.resolution = 8; % bits

%% M-PAM
pulse_shape = select_pulse_shape('rect', sim.ros.txDSP);
% pulse_shape = select_pulse_shape('rrc', sim.ros.txDSP, 0.25, 8);
mpam = PAM(sim.M, sim.Rb, 'equally-spaced', pulse_shape);

% estimate number of symbols necessary to minimize padding
sim.Nsymb  = Npoints/sim.ros.txDSP;

%% Summary
rows = {'DAC sampling rate'; 'Number of samples'; 'Number of symbols'; 'Trigger frequency';...
    'Symbol rate'; 'PAM order'; 'Bit rate';
    'Oversampling ratio of transmitter DSP'};
Variables = {'Tx.DAC.fs'; 'Npoints'; 'Nsymb'; 'ftrigger'; ...
    'mpam.Rs'; 'mpam.M'; 'mpam.Rb';...
    'sim.ros.txDSP'};
Values = [Tx.DAC.fs/1e9; Npoints; sim.Nsymb; ftrigger/1e9;...
    mpam.Rs/1e9; mpam.M; mpam.Rb/1e9;...
    sim.ros.txDSP];
Units = {'GS/s'; ''; ''; 'GHz'; 'GBaud'; ''; 'Gbit/s'; ''};

table(Variables, Values, Units, 'RowNames', rows)

%% Generate data
dataTX = randi([0 mpam.M-1], [1 sim.Nsymb]);
% dataTX = S.dataTX;
dataTX(1:sim.Nzero) = 0;
dataTX(end-sim.Nzero+1:end) = 0;

%% Duobinary enconding
if isfield(sim, 'duobinary') && sim.duobinary
    mpamdb = mpam.set_levels(0:mpam.M-1, 0.5 + (0:mpam.M-2));    
    
    xd = mpamdb.mod(dataTX); % Modulated PAM signal
    
    xd = duobinary_encoding(xd);
    
    xd_enc = xd;
    
    predist = @(p) 2/pi*asin(sqrt(p));
    dist = @(v) sin(pi/2*v)^2;
    
    Vswing = 1.8;
    Pswing = dist(Vswing/2);
    DP = Pswing/(mpam.M-1);
    Pk = 0:DP:Pswing;
    Vkp = predist(Pk);
    Vk = [-Vkp(end:-1:2) Vkp];
    Vk = Vk/Vk(end);
        
    ximp = upsample(Vk(xd_enc+mpam.M), mpam.pulse_shape.sps);
    xk = filter(ones(1, mpam.pulse_shape.sps)/mpam.pulse_shape.sps, 1, ximp);
%     xk = xk - mean(xk);
else
    if isfield(sim, 'mzm_predistortion') && sim.mzm_predistortion
        Vswing = 0.35*2;
        Vbias = 0.47;
        mpamPredist = mpam.mzm_predistortion(Vswing, Vbias, true);
        mpamPredist = mpamPredist.unbias;
        xk = mpamPredist.signal(dataTX); % Generate signal at sim.ros.txDSP rate
    else
        mpam = mpam.unbias;
        xk = mpam.signal(dataTX); % Generate signal at sim.ros.txDSP rate
    end
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
% [p, q] = rat(Tx.DAC.fs/(mpam.Rs*sim.ros.txDSP));
% xkDAC = resample(xk, p, q);
% [~, t1] = freq_time(length(xk), mpam.Rs*sim.ros.txDSP);
% dt = 1/Tx.DAC.fs;
% t2 = 0:dt:t1(end);
% xkDAC = interp1(t1, xk, t2);
% if t2(end) > t1(end)
%     Next = (t2(end)-t1(end))/
%     fprintf('Sequence had to be extended by \n', length(xkDAC) - Npoints)
%     
%     
% elseif t2(end) < t1(end)
%     fprintf('xkDAC was padded with %d points\n', Npoints - length(xkDAC))
%     xkDAC = [xkDAC ;
% end

xkDAC = xk;

%% Generate trigger signal
[~, t] = freq_time(length(xkDAC), Tx.DAC.fs);
xtrig = sin(2*pi*ftrigger*t);
xmin = min(xtrig);
xmax = max(xtrig);
xamp = xmax - xmin;    
dx = xamp/(2^(Tx.DAC.resolution)-1);  
sig_range = xmin:dx:xmax;
partition = sig_range(1:end-1) + dx/2;
[~, xqtrig] = quantiz(xtrig, partition, 0:255); 

% if not(isInteger(log2(length(xkDAC))))
%     Np2 = ceil(log2(length(xkDAC)));
%     xkDAC = [xkDAC xkDAC(1:Np2-length(xkDAC))];
%     xqtrig  = [xtrig xtrig(1:Np2-length(xkDAC))];
%     warning('Cyclic extension of %d points was necessary in order to make DAC waveform power of 2', Np2-length(xkDAC))
% end

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
sim.Mct = Tx.DAC.ros*10; % oversampling ratio to emulate continuous time
sim.N = sim.Nsymb*sim.Mct; % number of points 

%% Time and frequency
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'
[sim.f, sim.t] = freq_time(sim.N, sim.fs);

% DAC filter
% Tx.DAC.filt = design_filter('bessel', 5, 15e9/(sim.fs/2)); % DAC analog frequency response
fit3 = [-0.0013 0.5846 1.5859];
Tx.DAC.filt.H = @(f) 10.^(-polyval([-0.0013 0.5846 1.5859], abs(sim.fs*f/1e9))/20)...
    /10.^(-polyval([-0.0013 0.5846 1.5859], 0)/20);
% Tx.DAC.filt.H = Tx.DAC.filt.H(sim.f/sim.fs)/interp1(sim.f, Tx.DAC.filt.H(sim.f/sim.fs), 0);
% Tx.Mod.BW = interp1(Tx.Mod.H(sim.f > 0), sim.f(sim.f > 0),  0.5);

xt = dac(xk, Tx.DAC, sim, true);

figure(738)
subplot(211), box on, hold on
Np = 100;
plot([xq(end-Np:end) xq(1:Np)])
a = axis;
plot((Np+2)*[1 1], a(3:4), ':k')
set(gca, 'xtick', Np+2)
set(gca, 'xticklabel', 't = 0')
ylabel('xq')
legend('Signal', 'Transition')
title('Transmitted signal continuity')
subplot(212), box on, hold on
plot([xqtrig(end-Np:end) xqtrig(1:Np)])
a = axis;
plot((Np+2)*[1 1], a(3:4), ':k')
ylabel('xqtrig')
set(gca, 'xtick', Np+2)
set(gca, 'xticklabel', 't = 0')
legend('Trigger signal', 'Transition')
title('Trigger signal continuity')

% Plot
if sim.duobinary
    Ntraces = 500;
    Ndiscard = 100;
    Nstart = Ndiscard*sim.Mct+1;
    Nend = min(Nstart + Ntraces*2*sim.Mct, length(xt));
    figure(302), clf
    eyediagram(abs(sin(Vswing*pi/2*xt(Nstart:Nend))).^2, 2*sim.Mct)
    title('DAC output eye diagram')
    drawnow
end

%% Save file
filename = sprintf('data/waveforms/pam%d_%s_Rb=%dGbps', mpam.M, mpam.pulse_shape.type, round(mpam.Rb/1e9));
if sim.preemph
    filename = [filename '_preemph'];
end

if sim.duobinary
    filename = [filename '_duobin'];
elseif sim.mzm_predistortion
    filename = [filename '_predist'];
end

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
    fsrfDataJess = JESS_RAMFill_General(fsrfDataJess, xq, xqtrig, zeros(size(xq)), zeros(size(xq)));
end