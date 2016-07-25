function BER = DQPSK_DSP_BER_qsub(fiberLengthKm, Modulator, ModBWGHz, EqNtaps, linewidthKHz, fOffsetGHz, ros, ENOB)
%% Simulation of DSP-based coherent detection system using DQPSK
% - fiberLength: fiber length in km
% - Modulator: either 'EOM' or 'SiPhotonics'
% - ModBW: modulator bandwidth in GHz. Only used when Modulator == 'SiPhotonics'
% - EqNtaps: number of taps of adaptive equalizer
% - linewidth: transmitter and LO laser linewidth in kHz
% - fOffset: frequency offset of the LO laser with respect to transmitter
% laser in GHz. 
% - ros: oversampling ratio of DSP. Must be passed as a fraction e.g., 5/4
% - ENOB: effective number of bits

addpath DSP/
addpath f/
addpath ../f/
addpath ../apd/
addpath ../soa/

ros = eval(ros);
filename = sprintf('results/DQPSK_DSP_BER_L=%skm_%s_BW=%sGHz_Ntaps=%staps_nu=%skHz_fOff=%sGHz_ros=%d_ENOB=%s',...
        fiberLengthKm, Modulator, ModBWGHz, EqNtaps, linewidthKHz, fOffsetGHz, round(100*ros), ENOB)

% convert inputs to double (on cluster inputs are passed as strings)
if ~all(isnumeric([fiberLengthKm ModBWGHz EqNtaps linewidthKHz fOffsetGHz ENOB]))
    fiberLength = 1e3*str2double(fiberLengthKm);
    ModBW = 1e9*str2double(ModBWGHz);
    EqNtaps = round(str2double(EqNtaps));
    linewidth = 1e3*str2double(linewidthKHz);
    fOffset = 1e9*str2double(fOffsetGHz);
    ENOB = round(str2double(ENOB));
end

%% Simulation launched power swipe
Tx.PlaunchdBm = -38:-28;
% Tx.PlaunchdBm = -25;

%% ======================== Simulation parameters =========================
sim.Nsymb = 2^15; % Number of symbols in montecarlo simulation
sim.Mct = 10;    % Oversampling ratio to simulate continuous time 
sim.ros.rxDSP = ros;
sim.BERtarget = 1.8e-4; 
sim.Ndiscard = 1024; % number of symbols to be discarded from the begining and end of the sequence 
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation
sim.Rb = 2*112e9; % Bit rate
sim.Npol = 2;                                                              % number of polarizations
sim.Modulator = 'SiPhotonics';                                             % Modulator bandwidth limitation: {'MZM': limited by loss and velocity mismatch, 'SiPhotonics' : limited by parasitics (2nd-order response)}
sim.pulse_shape = select_pulse_shape('rect', sim.Mct);                     % pulse shape
sim.ModFormat = DPSK(4, sim.Rb/sim.Npol, sim.pulse_shape);                 % M-DPSK modulation format

% Simulation control
sim.RIN = true; 
sim.PMD = false;
sim.phase_noise = (linewidth ~= 0);
sim.preAmp = false;
sim.quantiz = not(isinf(ENOB));
sim.stopWhenBERreaches0 = true;                                            % whether to stop simulation after counter BER reaches 0
sim.save = true;

%% Plots
Plots = containers.Map();                                                   % List of figures 
Plots('BER')                  = 1; 
Plots('Eye diagram') = 0;
Plots('Channel frequency response') = 0;
Plots('Constellations') = 0;
Plots('Diff group delay')       = 0;
Plots('Phase tracker') = 0;
Plots('EPLL phase error') = 0;
Plots('Time recovery') = 0;
Plots('Phase error variance') = 0;
sim.Plots = Plots;
sim.shouldPlot = @(x) sim.Plots.isKey(x) && sim.Plots(x);

%% Time and frequency
sim.Rs = sim.ModFormat.Rs; % symbol rate
sim.fs = sim.ModFormat.Rs*sim.Mct;  % sampling frequency for 'continuous-time' simulation
[sim.f, sim.t] = freq_time(sim.N, sim.fs);

%% ===================== Transmitter Electric Filter ====================== 
Tx.filt = design_filter('bessel', 5, 0.7*sim.ModFormat.Rs/(sim.fs/2));      % design_filter(type, order, normalized cutoff frequency)
Tx.Delx  = 0;                                                               % Delay of x pol. in pol. mux. (symbol intervals)
Tx.Dely  = 0;                                                               % Delay of y pol. in pol. mux. (symbol intervals)

%% ========================== Transmitter Laser =========================== 
% Laser constructor: laser(lambda (nm), PdBm (dBm), RIN (dB/Hz), linewidth (Hz), frequency offset (Hz))
% lambda : wavelength (nm)
% PdBm : output power (dBm)
% RIN : relative intensity noise (dB/Hz)
% linewidth : laser linewidth (Hz)
% freqOffset : frequency offset with respect to wavelength (Hz)
Tx.Laser = laser(1250e-9, 0, -150, linewidth, 0);

%% ============================= Modulator ================================
if strcmpi(sim.Modulator, 'MZM') 
    %% Mach-Zehnder (limited by velocity mismatch and loss)
    % Mod = mzm_frequency_response(ratio : velocity mismatch, L : iteractive length in m, f: frequency vector, verbose: whether to plot results)
    Tx.Mod = mzm_frequency_response(0.98, 0.05, sim.f, true);
elseif strcmpi(sim.Modulator, 'SiPhotonics') 
    %% Si Photonics (limited by parasitics, 2nd-order response)
    Tx.Mod.BW = ModBW;
    Tx.Mod.fc = Tx.Mod.BW/sqrt(sqrt(2)-1); % converts to relaxation frequency
    Tx.Mod.grpdelay = 2/(2*pi*Tx.Mod.fc);  % group delay of second-order filter in seconds
    Tx.Mod.H = exp(1j*2*pi*sim.f*Tx.Mod.grpdelay)./(1 + 2*1j*sim.f/Tx.Mod.fc - (sim.f/Tx.Mod.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
end

%% =============================== Fiber ==================================
% Constructor: fiber(L, att(lamb) (optional), D(lamb) (optional)) 
% L : fiber length (m)
% att(lamb) : function handle of attenuation (att) at wavelength (lamb),
% deafault is att(lamb) = 0 dB/km
% D(lamb) : function handle of dispersion (D) at wavelength (lamb) in ps/(kmnm),
% default is D(lamb) = SSMF with lamb0 @ 1310 ps/(kmnm)
Fiber = fiber(fiberLength);
Fiber.PMD = sim.PMD;                                                       % whether to similate PMD
Fiber.meanDGDps = 0.1;                                                     % Mean DGD (ps)
Fiber.PMD_section_length = 1e3;                                            % Controls number of sections to simulate PMD (m)

%% ======================== Optical Amplifier =============================
% Constructor: soa(GaindB, NF, lamb, maxGain (optional))
% GaindB : Gain in dB
% NF : Noise figure in dB
% lamb : wavelength in m
% maxGain = maximum amplifier gain in dB, default is Inf
Amp = soa(20, 7, Tx.Laser.lambda);
% Note: class soa can be used for any amplifier, since it essentially 
% characterizes the amplifier in terms of gain and noise figure only

%% ======================= Local Oscilator ================================
Rx.LO = Tx.Laser;                                                          % Copy parameters from TX laser
Rx.LO.PdBm = 15;                                                           % Total local oscillator power (dBm)
Rx.LO.freqOffset = fOffset;                                                    % Frequency shift with respect to transmitter laser in Hz

%% ============================ Hybrid ====================================
% polarization splitting --------------------------------------------------
Rx.PolSplit.sig  = 'PBS';                                                   % pbs: polarization beamsplitter
Rx.PolSplit.LO   = 'PBS';                                                   % 3dB: 3-dB coupler     
Rx.PolSplit.Rext = 30;                                                      % PBS extinction ratio (dB), default = 30
Rx.PolSplit.R3dB = 1/2;                                                     % power splitting ratio of nominally 3-dB coupler (system performance should be insensitive to this parameter)
% 90-degree hybrid, same parameter for two polarizations ------------------
Rx.Hybrid.fS = 0.5;                                                         % power splitting ratio for signal coupler (W/W), default = 0.5
Rx.Hybrid.fL = 0.5;                                                         % power splitting ratio for LO coupler (W/W), default = 0.5
Rx.Hybrid.fI = 0.5;                                                         % power splitting ratio for in-phase coupler (W/W), default = 0.5
Rx.Hybrid.fQ = 0.5;                                                         % power splitting ratio for quadrature coupler (W/W), default = 0.5
Rx.Hybrid.tauIps = 0;                                                       % delay in in-phase branch (ps), default = 0
Rx.Hybrid.tauQps = 0;                                                       % delay in quadrature branch (ps), default = 0
Rx.Hybrid.phiI01deg = 0;                                                   % d.c. phase shift in I branch of pol. 1 (degrees), default = 0
Rx.Hybrid.phiQ01deg = 0;                                                    % d.c. phase shift in Q branch of pol. 1 (degrees), default = 0
Rx.Hybrid.phiI02deg = 0;                                                   % d.c. phase shift in I branch of pol. 2 (degrees), default = 0
Rx.Hybrid.phiQ02deg = 0;                                                    % d.c. phase shift in Q branch of pol. 2 (degrees), default = 0

%% ============================= Photodiodes ==============================
% Constructor: pin(R, Id, BW (optional))
% R : responsivity in A/W
% Id : dark current in A
% BW : bandwidth, default = Inf. Frequency response is a
% first-order filter with bandwidth BW.
Rx.PD = pin(1, 10e-9);

%% ======================== Transimpedance Amplifier ======================
Rx.N0 = (30e-12)^2;                                                        % One-sided thermal noise PSD per real dimension
% Note: ADC filter includes all the frequecy responses from the receiver

%% ================================= ADC ==================================
Rx.ADC.ENOB = ENOB;                                                           % effective number of bits
Rx.ADC.rclip = 0;                                                          % clipping ratio: clipped intervals: [xmin, xmin + xamp*rclip) and (xmax - xamp*rclip, xmax]
Rx.ADC.ros = sim.ros.rxDSP;                                                      % oversampling ratio with respect to symbol rate 
Rx.ADC.fs = Rx.ADC.ros*sim.Rs;                                             % ADC sampling rate
Rx.ADC.filt = design_filter('butter', 5, 0.5*Rx.ADC.fs/(sim.fs/2));            % design_filter(type, order, normalized cutoff frequency)
% ADC filter should include all filtering at the receiver: TIA,
% antialiasing, etc.

%% ================================= DSP ==================================
%% Equalization transmitter filter, CD, PMD, etc
Rx.AdEq.type  = 'CMA';                                                     % Adaptive equalization type: 'LMS' or 'CMA'. LMS only works when there's no phase noise or frequency offset
Rx.AdEq.structure = '2 filters';                                           % Structure: '2 filters' or '4 filters'. If '4 filters' corresponds to tranditional implementation, '2 filters' is simplified for short reach
Rx.AdEq.Ntrain = 1e4;                                                    % Number of symbols used for training (only used if LMS)
Rx.AdEq.mu = 1e-3;                                                         % Adaptation rate 
Rx.AdEq.Ntaps = EqNtaps;                                                         % Number of taps for each filter
Rx.AdEq.ros = sim.ros.rxDSP;                                                     % Oversampling ratio

%% Frequency recovery
% Used in DQPSK receiver and QPSK with feedforward carrier recovery
% Frequency recovery algorithm uses 4th power method to estimate frequency
% offset. Hence, it only works for quaternary modulation formats
Rx.FreqRec.mu = 5e-4;                                                % Adaptation rate. A vector means that gear shifting is used
Rx.FreqRec.Ntrain = 1e4;

%% Generate summary
coherent_simulation_summary(sim, Tx, Fiber, Rx);

%% Runs simulation
BER = ber_coherent_dsp(Tx, Fiber, Rx, sim);

%% Save results
if sim.save   
    % delete large variables
    sim = rmfield(sim, 'f');
    sim = rmfield(sim, 't');
    Tx.Mod = rmfield(Tx.Mod, 'H');    
    save(filename)
end
