function [BER, SNRdB] = QPSK_BER_qsub(fiberLength, Modulator, ModBW, CPRAlgorithm, CPRtaps, linewidth, fOffset, ros, ENOB)
%% Estimate BER of a QPSK system for parameters
% - fiberLength: fiber length in km
% - Modulator: either 'MZM' or 'SiPhotonics'
% - ModBW: modulator bandwidth in GHz. Only used when Modulator == 'SiPhotonics'
% - CPRAlgorithm: carrier phase recovery (CPR) algorithm. Either 'DPLL' or
% 'feedforward'
% - CPRtaps: Number of taps in CPR algorithm. Only used if CPRAlgorithm == 'feedforward'
% - linewidth: transmitter and LO laser linewidth in kHz
% - fOffset: frequency offset of the LO laser with respect to transmitter
% laser in GHz. 
% - ros: oversampling ratio of DSP
% - ENOB: effective number of bits

addpath DSP/
addpath f/
addpath ../f/
addpath ../apd/
addpath ../soa/
addpath ../mpam/

ros = eval(ros);
filename = sprintf('results/QPSK_BER_L=%skm_%s_BW=%sGHz_%s_%staps_nu=%skHz_fOff=%sGHz_ros=%d_ENOB=%s',...
        fiberLength, Modulator, ModBW, CPRAlgorithm, CPRtaps, linewidth, fOffset, round(100*ros), ENOB);

% convert inputs to double (on cluster inputs are passed as strings)
if ~all(isnumeric([fiberLength ModBW CPRtaps linewidth fOffset ENOB]))
    fiberLength = 1e3*str2double(fiberLength);
    ModBW = 1e9*str2double(ModBW);
    CPRtaps = round(str2double(CPRtaps));
    linewidth = 1e3*str2double(linewidth);
    fOffset = 1e9*str2double(fOffset);
    ENOB = round(str2double(ENOB));
end

%% ======================== Simulation parameters =========================
sim.Nsymb = 2^15;                                                          % Number of symbols in montecarlo simulation
sim.ros = ros;                                                             % Oversampling ratio of DSP
sim.Mct = 8*sim.ros;                                                       % Oversampling ratio to simulate continuous time 
sim.BERtarget = 1.8e-4;                                                    % Target BER
sim.Ndiscard = 1e3;                                                        % number of symbols to be discarded from the begining and end of the sequence 
sim.N = sim.Mct*sim.Nsymb;                                                 % number points in 'continuous-time' simulation
sim.Rb = 2*112e9;                                                          % Bit rate
sim.ModFormat = 'QAM';                                                     % either 'QAM' or 'DPSK'
sim.M    = 4;                                                              % QAM order
sim.pulse = 'NRZ';                                                         % Transmitter pulse shape ('NRZ')
sim.Npol = 2;                                                              % number of polarizations (currently ignored in generating the signal)
sim.Modulator = Modulator;                                                 % Modulator type: {'EOM': limited by loss and velocity mismatch, 'SiPhotonics' : limited by parasitics (2nd-order response)}
sim.Rs = sim.Rb/(sim.Npol*log2(sim.M));                                    % Symbol Rate
sim.save = true;

% Simulation control
sim.RIN = true; 
sim.PMD = true;
sim.phase_noise = (linewidth ~= 0);
sim.preAmp = false;  % currently ignored
sim.quantiz = ~isinf(ENOB);

%% Time and frequency
sim.fs = sim.Rs*sim.Mct;  % sampling frequency in 'continuous-time'

dt = 1/sim.fs;
t = 0:dt:(sim.N-1)*dt;
df = 1/(dt*sim.N);
f = -sim.fs/2:df:sim.fs/2-df;

sim.t = t;
sim.f = f;

%% Plots
Plots = containers.Map();                                                   % List of figures 
Plots('BER')                  = 1; 
Plots('Equalizer')            = 0;
Plots('ChannelFrequencyResponse') = 0;
Plots('CarrierPhaseNoise')    = 0;
Plots('DiffGroupDelay')       = 0;
Plots('PhaseTracker')         = 0;
Plots('FrequencyEstimation')  = 0; 
sim.Plots = Plots;

%% ===================== Transmitter Electric Filter ====================== 
Tx.filt = design_filter('bessel', 5, 0.7*sim.Rs/(sim.fs/2));               % design_filter(type, order, normalized cutoff frequency)
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
    % K. P. Ho, Phase-Modulated Optical Communication Systems. New York: Springer, 2005.
    % Parameters are the taken from Milad's code to result in a MZM with
    % bandwidth of about 50 GHz
    omega = 2*pi*f;
    c = 3e8;
    Mod.Vpi = 4.5;                                                         % Modulator switching voltage (V)
    ratio      = 0.98;                                                     % velocity mismatch ratio         
    n_r        = 2.15;                                                     % refractive index of coplanar-waveguide (CPW) for TM input light
    n_m        = n_r*ratio;                                                % 
    Mod.d_12   = (n_m -n_r)/c;                                             % Velocity mismatch difference between the optical and electrical waveguides
    Mod.a      = 0.01*100*sqrt(abs(omega)/2/pi*1e-9);                      % Microwave attenuation coefficient (electrode loss)
    Mod.L      = 0.05;                                                     % Interaction length 
    Mod.Hel    = (1-exp(Mod.L*(-Mod.a + 1j*Mod.d_12*omega)))./...          % Freq. response of the optical modulator limited by elec. loss and velocity mismatch
                    (Mod.a-1j*Mod.d_12*omega)/Mod.L;
    Mod.Hel(isnan(Mod.Hel)) = 1;                                           % Mod.Hel(f=0) is NaN  
elseif strcmpi(sim.Modulator, 'SiPhotonics') 
    %% Si Photonics (limited by parasitics, 2nd-order response)
    Mod.BW = ModBW;
    Mod.fc = Mod.BW/sqrt(sqrt(2)-1); % converts to relaxation frequency
    Mod.grpdelay = 2/(2*pi*Mod.fc);  % group delay of second-order filter in seconds
    Mod.Hel = exp(1j*2*pi*f*Mod.grpdelay)./(1 + 2*1j*f/Mod.fc - (f/Mod.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
end
Tx.Mod = Mod;                                                              % optical modulator

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
Rx.LO.PdBm = 9.8;                                                           % Total local oscillator power (dBm)
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
Rx.Hybrid.phiI01deg = 90;                                                   % d.c. phase shift in I branch of pol. 1 (degrees), default = 90
Rx.Hybrid.phiQ01deg = 0;                                                    % d.c. phase shift in Q branch of pol. 1 (degrees), default = 0
Rx.Hybrid.phiI02deg = 90;                                                   % d.c. phase shift in I branch of pol. 2 (degrees), default = 90
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
Rx.ADC.ros = sim.ros;                                                      % oversampling ratio with respect to symbol rate 
Rx.ADC.fs = Rx.ADC.ros*sim.Rs;                                             % ADC sampling rate
Rx.ADC.filt = design_filter('cheby1', 5, 0.5*Rx.ADC.fs/(sim.fs/2));            % design_filter(type, order, normalized cutoff frequency)
% ADC filter should include all filtering at the receiver: TIA,
% antialiasing, etc.

%% ================================= DSP ==================================
%% Equalization transmitter filter, CD, PMD, etc
Rx.AdEq.type  = 'CMA';                                                     % Adaptive equalization type: 'LMS' or 'CMA'. LMS only works when there's no phase noise or frequency offset
Rx.AdEq.structure = '2 filters';                                           % Structure: '2 filters' or '4 filters'. If '4 filters' corresponds to tranditional implementation, '2 filters' is simplified for short reach
Rx.AdEq.Ntrain = 1e4;                                                    % Number of symbols used for training (only used if LMS)
Rx.AdEq.mu = 1e-3;                                                         % Adaptation rate 
Rx.AdEq.Ntaps = 3;                                                         % Number of taps for each filter
Rx.AdEq.ros = sim.ros;                                                     % Oversampling ratio

%% Carrrier phase recovery
% Only used if sim.ModFormat = 'QAM'
Rx.CPR.type = CPRAlgorithm;                                                % Carrier phase recovery: 'DPLL' (digital phase-locked loop) or 'feedforward'
Rx.CPR.phaseEstimation = 'NDA';                                             % Phase estimation method: 'dd' = decision-directed for either DPLL or feedfoward; 'nd' = non-data-aided (only for feedforward); '4th power' (only for DPLL)
Rx.CPR.Delay = 0;                                                       % Delay in number of symbols due to pipelining and parallelization
Rx.CPR.Ntrain = 1e4;                                                     % Number of symbols used for training. This training starts when equalization training is done
% Carrier phase recovery parameters for 'feedforward'
Rx.CPR.FilterType = 'FIR';                                                 % Filter type: 'FIR' or 'IIR'
Rx.CPR.structure = '2 filters';                                             % structure of feedforward employing DD and FIR filter: {'1 filter', '2 filter'}
Rx.CPR.Ntaps = CPRtaps;                                                          % Maximum number of taps of filter 
Rx.CPR.NDAorder = 'reverse';                                               % Order of phase estimation (PE) and filtering in NDA FF:{'direct': PE followed by filtering, 'reverse': filtering followed by PE}
% Carrier phase recovery parameters for 'DPLL'
Rx.CPR.CT2DT = 'bilinear';                                                 % method for converting continuous-time loop filter to discrete time: {'bilinear', 'impinvar'}
Rx.CPR.csi = 1/sqrt(2);                                                    % damping coefficient of second-order loop filter
Rx.CPR.wn = 2*pi*0.8e9;                                                    % relaxation frequency of second-order loop filter: optimized using optimize_PLL.m
% Phase tracking stage
% This only compensates for constant phase offsets, which may be required
% since feedforward or DPLL may converge to a rotated constellation when
% the frequency offset is large
% Phase tracking is adpated using LMS
Rx.PT.mu = 0.001;
Rx.PT.Ntrain = 0.5e4;

%% Frequency recovery
% Only used if sim.ModFormat = 'DPSK', since for 'QAM' feedforward or DPLL
% can compensate frequency offset as well
% Frequency recovery algorithm uses 4th power method to estimate frequency
% offset. Hence, it only works for DQPSK
Rx.FreqRec.mu = [0.0005 0];                                                 % Adaptation rate. A vector means that gear shifting is used
Rx.FreqRec.muShift = 2e4;                                  % Controls when gears are shifted. From 1:muShift(1) use mu(1), from muShift(1)+1:muShift(2) use mu(2), etc
Rx.FreqRec.Ntrain = 2e4;

Tx.PlaunchdBm = -35:0.5:-29;

% [Nsum, Nmult] = calcDSPOperations(Rx, sim)

[BER, SNRdB] = ber_coherent(Tx, Fiber, Rx, sim)

if sim.save   
    % delete large variables
    clear t f Mod omega Fiber c
    sim = rmfield(sim, 'f');
    sim = rmfield(sim, 't');
    Tx.Mod = rmfield(Tx.Mod, 'Hel');
    if isfield(Tx.Mod, 'a')
        Tx.Mod = rmfield(Tx.Mod, 'a');
    end
    
    save(filename)
end
