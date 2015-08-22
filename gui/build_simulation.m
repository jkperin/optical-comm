function [mpam, ofdm1, tx, fiber1, soa1, apd1, rx, sim] = build_simulation(h)
%% Extract user-inputted values from GUI and build required classes and structs

% Auxiliary functions
getString = @(h) get(h, 'String');
getValue = @(h) str2double(get(h, 'String'));
getLogicalValue = @(h) logical(get(h, 'Value'));

%% Simulation parameters
sim.Nsymb = eval(getString(h.Nsymb)); % Number of symbols in montecarlo simulation
sim.Mct = getValue(h.Mct);     % Oversampling ratio to simulate continuous time (must be odd) 
if mod(sim.Mct, 2) == 0
    warning('Oversampling ratio for continuous time must be odd. Mct = %d provided, using Mct = %d instead', sim.Mct, sim.Mct+1);
    sim.Mct = sim.Mct+1;
end

sim.BERtarget = eval(getString(h.BER)); 

sim.shot = getLogicalValue(h.check.shot); % include shot noise in montecarlo simulation (always included for pin and apd case)
sim.RIN = getLogicalValue(h.check.rin); % include RIN noise in montecarlo simulation
sim.quantiz = getLogicalValue(h.check.quantiz);
sim.ENOB = getValue(h.ENOB);

sim.verbose = false; % show stuff
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation
sim.L = getValue(h.Lseq);  % de Bruijin sub-sequence length (ISI symbol length)

%% Modulation
modulation = getOption(h.popup.modulation);
switch modulation
    case 'M-PAM'
        switch getOption(h.level)
            case 'Equal'
                level_spacing = 'equally-spaced'; % M-PAM level spacing: 'uniform' or 'non-uniform'
            case 'Optimized'
                level_spacing = 'optimized';
            otherwise
                error('Level spacing option not yet implemented')
        end
        M = getValue(h.M);
        Rb = 1e9*getValue(h.Rb);
        
        
        str = get(h.pshape, 'String');
        if strcmp(str{get(h.pshape, 'Value')}, 'Rectangular')
            pshape = @(n) double(n >= 0 & n < sim.Mct); % pulse shape
        else
            error('Root-Raised Cosine not yet implemented')
        end
        
        mpam = PAM(M, Rb, level_spacing, pshape);
        
        sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'
        
    case 'DMT/OFDM'
        sim.type = getOption(h.palloc);
        switch sim.type
            case 'Preemphasis'
                sim.type = 'preemphasis';
            case 'Opt. Bit Loading & Power allocation'
                sim.type = 'palloc';
        end
        
        sim.rclip = 10^(getValue(h.rclip)/20);
        
        Nc = getValue(h.Nc);
        Nu = getValue(h.Nu);
        CS = getValue(h.M);
        Rb = 1e9*getValue(h.Rb);
        
        ofdm1 = ofdm(Nc, Nu, CS, Rb);
        
    case 'M-CAP'
        error('M-CAP not yet implemented')
    otherwise
        error('Invalid Option')
end

%% Transmitter
tx.PtxdBm = eval(get(h.Ptx, 'String'));
tx.lamb = 1e-9*eval(getString(h.lamb));
tx.kappa = 1;

if getLogicalValue(h.check.chirp)
    tx.alpha = getValue(h.chirp);
end

if sim.RIN
    tx.RIN = getValue(h.rin);   % RIN in dB/Hz. Only used if sim.RIN is true
    if getLogicalValue(h.check.rin_shape)
        g = getValue(h.rin_shape);
        tx.RIN_variation = h.default.RIN_variation;
        tx.RIN_bw = h.default.RIN_bw;
        tx.RIN_shape = @(f, fr) (f.^2 + (g/(2*pi)).^2)./((fr^2-f.^2).^2 + f.^2*(g/(2*pi)).^2);
    end
else
    tx.RIN = -Inf;
end

if getLogicalValue(h.check.rex)
    tx.rexdB = -abs(getValue(h.rex));   % exctinction ratio in dB. Defined as Pmin/Pmax
else
    tx.rexdB = -Inf;
end

% Modulator frequency response
tx.modulator.Fc = 1e9*eval(get(h.fc, 'String')); % modulator cut off frequency

%% Fiber
L = 1e3*eval(getString(h.L));
att = getValue(h.att);
D = 1e-6*getValue(h.D);

if length(L) ~= 1
    sim.fiberL = L;
end

fiber1 = fiber(L(1), @(l) att, @(l) D);

%% Receiver
rx.R = getValue(h.R); % responsivity
rx.NEP = getValue(h.NEP)*1e-12;
rx.N0 = (rx.NEP).^2; % thermal noise psd
rx.Sth = rx.N0/2;
rx.Id = getValue(h.Id)*1e-9; % dark current

% Electric Lowpass Filter
switch getOption(h.rxfilterType)
    case 'Matched'
        rx.filter.type = 'matched';
    case 'Gaussian'
        rx.filter.type = 'gaussian';
    case 'Bessel'
        rx.filter.type = 'bessel';
    case 'Butterworth'
        rx.filter.type = 'butter';
    otherwise
        error('Invalid receiever filter option')
end

rx.filterN = getValue(h.rxfilterN);
rx.filterBw = 1e9*getValue(h.rxfilterBw);
   
%% Equalization
str = get(h.eq_type, 'String');
val = get(h.eq_type, 'Value');
rx.eq.type = str{val};
rx.eq.ros = getValue(h.eq_ros);
rx.eq.Ntaps = getValue(h.eq_taps);

if strfind(rx.eq.type, 'Adaptive')
    rx.eq.adaptive = true;
    rx.eq.mu = getValue(h.eq_mu);
    rx.eq.Ntrain = getValue(h.eq_Ntrain);
else
    rx.eq.adaptive = false;
end

%% ADS co-simulation
sim.ads.cosim = getLogicalValue(h.check.ads);
sim.ads.model = getOption(h.ads);
sim.ads.eyediagram = getLogicalValue(h.check.ads_eye);
sim.ads.risetime = str2double(strrep(getOption(h.trise), ' ps', ''));

switch sim.ads.model
    case 'DFB 25C'
        sim.ads.Temp = 25;
        sim.ads.Plevels = [3.5 6.2 9 11.8]; % rough estimates of the levels
        sim.ads.Pthresh = [5 7.8 10.6]; % rough estimates of the decision thresholds
    case 'DFB -5C'
        sim.ads.Temp = -5;
        sim.ads.Plevels = [3 6 9 12.5]; % rough estimates of the levels
        sim.ads.Pthresh = [4.3 7.6 11]; % rough estimates of the decision thresholds
    case 'DFB 75C'
        sim.ads.Temp = 75;
        sim.ads.Plevels = [4 6.8 9.8 12.2]; % rough estimates of the levels
        sim.ads.Pthresh = [5.35 8.25 11.3]; % rough estimates of the decision thresholds
    otherwise 
        error('Invalid ADS model.')
end

sim.ads.filename = sprintf('ADS_DFB_%dC_%dps.mat', sim.ads.Temp, sim.ads.risetime);

%% System simulation
switch getOption(h.popup.system)
    case 'Basic'
        apd1 = [];
        soa1 = [];
        if strcmp(modulation, 'M-PAM')
            ofdm1 = [];
        elseif strcmp(modulation, 'DMT/OFDM')
            mpam = [];
        end
    case 'SOA'     
        if strcmp(getOption(h.optfiltType), 'Fabry-Perot')
            filterType = 'lorentzian';
        else
            filterType = 'fbg';
        end
        
        Bopt = 1e9*getValue(h.Bopt);
        sim.M = ceil(Bopt/mpam.Rs); % Ratio of optical filter BW and electric filter BW (must be integer)
        sim.Me = max(2*sim.M, 16);  % number of used eigenvalues
        sim.polarizer = getLogicalValue(h.check.polarizer);
        
        % Optical Bandpass Filter
        rx.optfilt = design_filter(filterType, 0, Bopt/(sim.fs/2));
        
        GsoadB = getValue(h.Gsoa);
        Fn = getValue(h.Fn);
        
        soa1 = soa(GsoadB, Fn, tx.lamb, Inf); 
                       
        apd1 = [];
        
        ofdm1 = [];
                     
    case 'APD'
        if getLogicalValue(h.check.GBw)
            GBw = 1e9*getValue(h.GBw);
        else
            GBw = Inf;
        end
        
        ka = getValue(h.ka);
        GapddB =  getValue(h.Gapd);
        
        apd1 = apd(GapddB, ka, GBw, rx.R, rx.Id);
        
        sim.OptimizeGain = ~getLogicalValue(h.check.Gapd);
        
        soa1 = [];
        
        ofdm1 = [];
end
   
end

function str = getOption(h)
list = get(h, 'String');
str = list{get(h, 'Value')};
end