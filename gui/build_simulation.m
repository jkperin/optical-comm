function [mpam, ofdm1, tx, fiber1, soa1, apd1, rx, sim] = build_simulation(h)

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
        if getOption(h.level)
            mpam.level_spacing = 'uniform'; % M-PAM level spacing: 'uniform' or 'non-uniform'
        else
            mpam.level_spacing = 'nonuniform';
        end
        mpam.M = getValue(h.M);
        mpam.Rb = 1e9*getValue(h.Rb);
        mpam.Rs = mpam.Rb/log2(mpam.M);
        
        str = get(h.pshape, 'String');
        if strcmp(str{get(h.pshape, 'Value')}, 'Rectangular')
            mpam.pshape = @(n) double(n >= 0 & n < sim.Mct); % pulse shape
        else
            error('Root-Raised Cosine not yet implemented')
        end
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
tx.lamb = 1e-9*getValue(h.lamb);
tx.kappa = 10^(getValue(h.kappa)/10);

if getLogicalValue(h.check.chirp)
    tx.alpha = getValue(h.chirp);
end

if sim.RIN
    tx.RIN = getValue(h.rin);   % RIN in dB/Hz. Only used if sim.RIN is true
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
att = @(lamb) getValue(h.att);
D = @(lamb) getValue(h.D);

if length(L) ~= 1
    sim.fiberL = L;
end

fiber1 = fiber(L(1), att, D);

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
    otherwise
        error('Invalid receiever filter option')
end

rx.filterN = getValue(h.rxfilterN);
rx.filterBw = 1e9*getValue(h.rxfilterBw);
   
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
        
        % Optical Bandpass Filter
        rx.optfilt = design_filter(filterType, 0, Bopt/(sim.fs/2));
        
        GsoadB = getValue(h.Gsoa);
        Fn = getValue(h.Fn);
        maxGsoadB = getValue(h.maxGsoa);
        
        soa1 = soa(GsoadB, Fn, tx.lamb, maxGsoadB); 
        
        sim.OptimizeGain = getLogicalValue(h.check.OptGsoa);
               
        apd1 = [];
        
        ofdm1 = [];
                     
    case 'APD'
        GBw = 1e9*getValue(h.GBw);
        ka = getValue(h.ka);
        GapddB =  getValue(h.Gapd);
        
        apd1 = apd(GapddB, ka, GBw, rx.R, rx.Id);
        
        sim.OptimizeGain = getLogicalValue(h.check.OptGapd);
        
        soa1 = [];
        
        ofdm1 = [];
end
   



end

function str = getOption(h)
list = get(h, 'String');
str = list{get(h, 'Value')};
end