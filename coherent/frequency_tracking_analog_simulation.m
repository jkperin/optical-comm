%% Simulation of analog coherent detection system
clear, clc

addpath f/
addpath analog/
addpath ../f/
addpath ../apd/
addpath ../soa/

%% ======================== Simulation parameters =========================
sim.Nsymb = 2^13; % Number of symbols in montecarlo simulation
sim.Mct = 10;    % Oversampling ratio to simulate continuous time 
sim.BERtarget = 1.8e-4; 
sim.Ndiscard = 256; % number of symbols to be discarded from the begining and end of the sequence 
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation
sim.Rb = 2*112e9; % Bit rate
sim.Npol = 2;                                                              % number of polarizations
sim.Modulator = 'SiPhotonics';                                             % Modulator bandwidth limitation: {'MZM': limited by loss and velocity mismatch, 'SiPhotonics' : limited by parasitics (2nd-order response)}
sim.pulse_shape = select_pulse_shape('rect', sim.Mct);                     % pulse shape
sim.ModFormat = QAM(4, sim.Rb/sim.Npol, sim.pulse_shape);                  % M-QAM modulation format

% Simulation control
sim.RIN = true; 
sim.PMD = false;
sim.phase_noise = ~true;
sim.preAmp = false;
sim.stopWhenBERreaches0 = true;                                            % whether to stop simulation after counter BER reaches 0

%% Plots
Plots = containers.Map();                                                   % List of figures 
Plots('Channel frequency response') = 0;
Plots('Constellations') = 1;
Plots('Time recovery') = 0;
Plots('Phase error variance') = 0;
Plots('Symbol errors') = 1;
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
Tx.Laser = laser(1250e-9, 0, -150, 200e3, 0);
Tx.PlaunchdBm = -20;

%% ============================= Modulator ================================
if strcmpi(sim.Modulator, 'MZM') 
    %% Mach-Zehnder (limited by velocity mismatch and loss)
    % Mod = mzm_frequency_response(ratio : velocity mismatch, L : iteractive length in m, f: frequency vector, verbose: whether to plot results)
    Tx.Mod = mzm_frequency_response(0.98, 0.05, sim.f, true);
elseif strcmpi(sim.Modulator, 'SiPhotonics') 
    %% Si Photonics (limited by parasitics, 2nd-order response)
    Tx.Mod.BW = 40e9;
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
Fiber = fiber(0e3);
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

% Frequency step
% Tstep = sim.N/2; % frequency step begins
% fstep = 0.5e9;
% Rx.LO.freqOffset = [zeros(1, Tstep-1) fstep*ones(1, sim.N-Tstep+1)];    % Frequency shift with respect to transmitter laser in Hz

% Frequency ramp
Tramp = sim.N/4; % frequency ramp begins
framp = 3e9; % frequency ramp slope in Hz/us
Rx.LO.freqOffset = [zeros(1, Tramp-1) 1e6*framp*(0:sim.N-Tramp)/sim.fs];    % Frequency shift with respect to transmitter laser in Hz     

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

%% ========================= Analog Components ============================
%% Receiver filter
Analog.filt = design_filter('butter', 5, 0.7*sim.ModFormat.Rs/(sim.fs/2));

%% Carrier phase recovery and components
% Carrier Phase recovery type: either 'OPLL' (not implemented), 'EPLL',
% and 'Feedforward'
Analog.CarrierPhaseRecovery = 'EPLL';
% CPRmethod: {'Costas': electric PLL based on Costas loop, which
% requires multiplications, 'logic': EPLL based on XOR operations, 
% '4th-power': based on raising signal to 4th power}
Analog.CPRmethod = 'logic';                                            

% If componentFilter is empty, simulations assume ideal devices
componentFilter = []; %design_filter('bessel', 1, 0.5*sim.Rs/(sim.fs/2));
componentRn = 60; % (Ohm) equivalent noise resistance obtained from 
% Huber, A. et al (2002). Noise model of InP-InGaAs SHBTs for RF circuit design. 
% IEEE Transactions on Microwave Theory and Techniques, 50(7), 1675–1682.
componentN0 = 4e-21*componentRn/pi;

% Adder
Analog.Adder.filt = componentFilter;
Analog.Adder.N0 = componentN0;

% Mixer
Analog.Mixer.filt = componentFilter;
Analog.Mixer.N0 = componentN0;

% ABS (full-wave rectifier)
Analog.ABS.filt = componentFilter;
Analog.ABS.N0 = componentN0;

% Logic
Analog.Logic.Vcc = 1;
Analog.Logic.N0 = componentN0;
Analog.Logic.filt = componentFilter;

% Comparator
Analog.Comparator.Vref = 0;
Analog.Comparator.Vcc = 1;
Analog.Comparator.N0 = componentN0;
Analog.Comparator.filt = componentFilter;

%% PLL loop filter parameters.
% Note: relaxation frequency is optimized at every iteration
Analog.Kdc = 1;                                                            % DC gain
Analog.csi = 1/sqrt(2);                                                    % damping coefficient of second-order loop filter
Analog.Delay = 0;                                                        % Additional loop delay in s (not including group delay from filters)

%% Feedforward additional components
Analog.Feedforward.freqDivDelay = 0;                                       % delay in samples
Analog.Feedforward.FreqDiv1.Mixer.filt = [];
Analog.Feedforward.FreqDiv1.Mixer.N0 = 0;
Analog.Feedforward.FreqDiv2.Mixer.filt = [];
Analog.Feedforward.FreqDiv2.Mixer.N0 = 0;

Analog.Feedforward.LPF.filt = design_filter('bessel', 5, 5e9/(sim.fs/2));
Analog.Feedforward.FreqDiv1.filt = design_filter('bessel', 5, 2e9/(sim.fs/2));
Analog.Feedforward.FreqDiv2.filt = design_filter('bessel', 5, 1e9/(sim.fs/2));

Rx.Analog = Analog;

%% ========================= Time recovery ================================
% Two types are supported: 'spectral-line'
% Spectral line method: nonlinearity (squarer) -> BPF or PLL
Rx.TimeRec.type = 'spectral-line-bpf';
Rx.TimeRec.type = 'spectral-line-pll';
Rx.TimeRec.type = 'none';
Rx.TimeRec.Mct = 30; % oversampling ratio of continuous time used in TimeRecovery
Rx.TimeRec.fs = Rx.TimeRec.Mct*sim.ModFormat.Rs;
Rx.TimeRec.squarerFilt = design_filter('fir', 32, sim.ModFormat.Rs/(Rx.TimeRec.fs/2)); 
% Note: use FIR filter to model squarerFilt so that group delay can be
% easily removed in order to preserve timing accuracy.
Rx.TimeRec.N0 = componentN0; % noise PSD of circuitry
Rx.TimeRec.Ndiscard = [2e3 512];

% Additional paramters for 'spectral-line-bpf'
BW = 1e9;
lpf = design_filter('butter', 5, BW/(Rx.TimeRec.fs/2));
[bpf.num, bpf.den] = iirlp2bp(lpf.num, lpf.den, BW/(Rx.TimeRec.fs/2), sim.ModFormat.Rs/(Rx.TimeRec.fs/2) + BW/(Rx.TimeRec.fs/2)*[-1 1]); % converts to BPF
Rx.TimeRec.bpf = bpf; % this will be filtered using filtfilt (zero-phase filtering)

% Additional paramters for 'spectral-line-pll'
Rx.TimeRec.csi = sqrt(2)/2; % damping
Rx.TimeRec.wn = 2*pi*1e9; % relaxation frequency of PLL
Rx.TimeRec.CT2DT = 'bilinear'; % continuous-time to discrete-time conversion method 
Rx.TimeRec.mixerFilt = []; % design_filter('butter', 5, (sim.ModFormat.Rs)/(Rx.TimeRec.fs/2)); % filtering operation performed before and after squaring
Rx.TimeRec.loopDelay = 10; % additional loop delay
    
%% Generate summary
coherent_simulation_summary(sim, Tx, Fiber, Rx);

%% Hold in range
frequency_tracking_coherent_analog_qpsk(Tx, Fiber, Rx, sim)
