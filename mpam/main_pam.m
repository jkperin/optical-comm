%% Calculate BER of IM-DD system using M-PAM
% - Equalization is done using a fractionally spaced linear equalizer
% Simulations include modulator, fiber, optical amplifier (optional) characterized 
% only by gain and noise figure, optical bandpass filter, antialiasing 
% filter, sampling, and linear equalization
clear, clc

addpath ../f % general functions
addpath ../apd % for PIN photodetectors

%% Transmit power swipe
Tx.PtxdBm = -30:-22; % transmitter power range
% Tx.PtxdBm = -13; % transmitter power range

%% Simulation parameters
sim.Rb = 56e9;    % bit rate in bits/sec
sim.Nsymb = 2^14; % Number of symbols in montecarlo simulation
sim.ros.txDSP = 2; % oversampling ratio transmitter DSP (must be integer). DAC samping rate is sim.ros.txDSP*mpam.Rs
% For DACless simulation must make Tx.dsp.ros = sim.Mct and DAC.resolution = Inf
sim.ros.rxDSP = 2; % oversampling ratio of receiver DSP. If equalization type is fixed time-domain equalizer, then ros = 1
sim.Mct = 2*6;      % Oversampling ratio to simulate continuous time. Must be integer multiple of sim.ros.txDSP and numerator of sim.ros.rxDSP
sim.BERtarget = 1e-4; 
sim.Ndiscard = 1024; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation
sim.Modulator = 'MZM'; % 'MZM' or 'DML'
 
%% Simulation control
sim.preAmp = true;
sim.preemphasis = false; % preemphasis to compensate for transmitter bandwidth limitation
sim.preemphRange = 25e9; % preemphasis range
sim.mzm_predistortion = 'levels'; % predistortion to compensate MZM nonlinearity {'none': no predistortion, 'levels': only PAM levels are predistorted, 'analog': analog waveform is predistorted (DEPRECATED)}
sim.RIN = true; % include RIN noise in montecarlo simulation
sim.phase_noise = true; % whether to simulate laser phase noise
sim.PMD = false; % whether to simulate PMD
sim.quantiz = ~true; % whether quantization is included
sim.stopSimWhenBERReaches0 = true; % stop simulation when counted BER reaches 0

% Control what should be plotted
sim.Plots = containers.Map();
sim.Plots('BER') = 1;
sim.Plots('DAC output') = 0;
sim.Plots('Optical eye diagram') = 1;
sim.Plots('Received signal eye diagram') = 0;
sim.Plots('Signal after equalization') = 1;
sim.Plots('Equalizer') = 1;
sim.Plots('Electronic predistortion') = 1;
sim.Plots('Adaptation MSE') = 0;
sim.Plots('Channel frequency response') = 0;
sim.Plots('OSNR') = 1;
sim.Plots('Received signal optical spectrum') = 0;
sim.Plots('PAM levels MZM predistortion') = 0;
sim.shouldPlot = @(x) sim.Plots.isKey(x) && sim.Plots(x);

%% Pulse shape
Tx.pulse_shape = select_pulse_shape('rect', sim.ros.txDSP);
% pulse_shape = select_pulse_shape('rrc', sim.ros.txDSP, 0.5, 6);

%% M-PAM
% PAM(M, bit rate, leve spacing : {'equally-spaced', 'optimized'}, pulse
% shape: struct containing properties of pulse shape 
mpam = PAM(4, sim.Rb, 'equally-spaced', Tx.pulse_shape);

%% Time and frequency
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'
[sim.f, sim.t] = freq_time(sim.N, sim.fs);

%% ==========================  Transmitter ================================
%% DAC
Tx.DAC.fs = sim.ros.txDSP*mpam.Rs; % DAC sampling rate
Tx.DAC.ros = sim.ros.txDSP; % oversampling ratio of transmitter DSP
Tx.DAC.resolution = 6; % DAC effective resolution in bits
Tx.DAC.rclip = 0;
Tx.DAC.filt = design_filter('bessel', 5, 30e9/(sim.fs/2)); % DAC analog frequency response
% juniperDACfit = [-0.0013 0.5846 0];
% Tx.DAC.filt.H = @(f) 1./(10.^(polyval(juniperDACfit, abs(f*sim.fs/1e9))/20)); % DAC analog frequency response

%% Laser
% Laser constructor: laser(lambda (nm), PdBm (dBm), RIN (dB/Hz), linewidth (Hz), frequency offset (Hz))
% lambda : wavelength (nm)
% PdBm : output power (dBm)
% RIN : relative intensity noise (dB/Hz)
% linewidth : laser linewidth (Hz)
% freqOffset : frequency offset with respect to wavelength (Hz)
Tx.Laser = laser(1550e-9, 0, -150, 0.2e6, 0);

%% Modulator
Tx.rexdB = -20;  % extinction ratio in dB. Defined as Pmin/Pmax
Tx.Mod.type = sim.Modulator; 
Tx.Vgain = 1; % Gain of driving signal
Tx.VbiasAdj = 1; % adjusts modulator bias
Tx.Mod.Vbias = 0.5; % bias voltage normalized by Vpi
Tx.Mod.Vswing = 1;  % normalized voltage swing. 1 means that modulator is driven at full scale
Tx.Mod.BW = 20e9; % DAC frequency response includes DAC + Driver + Modulator
Tx.Mod.filt = design_filter('bessel', 5, Tx.Mod.BW/(sim.fs/2));
Tx.Mod.H = Tx.Mod.filt.H(sim.f/sim.fs);
Tx.Mod.alpha = 0;

%% ============================= Fiber ====================================
% fiber(Length in m, anonymous function for attenuation versus wavelength
% (default: att(lamb) = 0 i.e., no attenuation), anonymous function for 
% dispersion versus wavelength (default: SMF28 with lamb0 = 1310nm, S0 = 0.092 s/(nm^2.km))
SMF = fiber(0e3); 
DCF = fiber(0, @(lamb) 0, @(lamb) -0.1*(lamb-1550e-9)*1e3 - 40e-6); 

Fibers = [SMF DCF];

linkAttdB = SMF.att(Tx.Laser.wavelength)*SMF.L/1e3...
    + DCF.att(Tx.Laser.wavelength)*DCF.L/1e3;

%% ========================== Amplifier ===================================
% Constructor: OpticalAmplifier(GaindB, NF, lamb, maxGain (optional))
% GaindB : Gain in dB
% NF : Noise figure in dB
% lamb : wavelength in m
% maxGain = maximum amplifier gain in dB, default is Inf
Rx.OptAmp = OpticalAmplifier(30, 5, Tx.Laser.lambda);
Rx.OptAmpOutPowerdBm = 0; % output power after amplifier
% Note: the amplifier here operates in the constant output power mode,
% where the output power after amplification is set to Rx.AmpOutPowerdBm

%% ============================ Receiver ==================================
%% Photodiodes
% pin(R: Responsivity (A/W), Id: Dark current (A), BW: bandwidth (Hz))
% PIN frequency response is modeled as a first-order system
Rx.PD = pin(1, 10e-9, Inf);

%% TIA-AGC
% One-sided thermal noise PSD
Rx.N0 = (30e-12).^2; 

%% Receiver DSP
Rx.filtering = 'antialiasing'; % {'antialiasing' or 'matched'}

%% ADC for direct detection case
Rx.ADC.fs = mpam.Rs*sim.ros.rxDSP;
Rx.ADC.ros = sim.ros.rxDSP;
Rx.ADC.filt = design_filter('butter', 5, (Rx.ADC.fs/2)/(sim.fs/2)); % Antialiasing filter
Rx.ADC.ENOB = 5; % effective number of bits. Quantization is only included if sim.quantiz = true and ENOB ~= Inf
Rx.ADC.rclip = 0;

%% Equalizer
% Terminology: TD = time domain, SR = symbol-rate, LE = linear equalizer
Rx.eq.ros = sim.ros.rxDSP;
Rx.eq.type = 'Adaptive TD-LE';
Rx.eq.Ntaps = 9; % number of taps
Rx.eq.mu = 1e-3; % adaptation ratio
Rx.eq.Ntrain = 4e3; % Number of symbols used in training (if Inf all symbols are used)
Rx.eq.Ndiscard = [5e3 1024]; % symbols to be discard from begining and end of sequence due to adaptation, filter length, etc

%% Generate summary
% generate_summary(mpam, Tx, Fibers, Rx.OptAmp, Rx, sim);

%% Run simulation
%% Calculate BER
Ptx = dBm2Watt(Tx.PtxdBm); % Transmitted power

ber.gauss = zeros(size(Ptx));
ber.count = zeros(size(Ptx));
OSNRdB = zeros(size(Ptx));
for k = 1:length(Ptx)
    Tx.Ptx = Ptx(k);
         
    % Montecarlo simulation
    [ber.count(k), ber.gauss(k), OSNRdB(k), Rx] = ber_pam_montecarlo(mpam, Tx, Fibers, Rx, sim);
    
    if sim.stopSimWhenBERReaches0 && ber.count(k) == 0
        break;
    end 
end

%% Plots
if sim.shouldPlot('BER') && length(ber.count) > 1
    figure(1), hold on, box on
    if sim.preAmp
        hline(1) = plot(OSNRdB, log10(ber.count), '-o', 'LineWidth', 2, 'DisplayName', 'Counted');
        hline(2) = plot(OSNRdB, log10(ber.gauss), '-', 'Color', get(hline(1), 'Color'), 'LineWidth', 2, 'DisplayName', 'Gaussian Approximation');
        hline(3) = plot(OSNRdB(OSNRdB ~= 0), log10(pam_ber_from_osnr(mpam.M, OSNRdB(OSNRdB ~= 0), Rx.noiseBW)), '--k', 'LineWidth', 2, 'DisplayName', 'Sig-spont & noise enhancement');
        hline(3) = plot(OSNRdB(OSNRdB ~= 0), log10(pam_ber_from_osnr(mpam.M, OSNRdB(OSNRdB ~= 0), mpam.Rs/2)), ':k', 'LineWidth', 2, 'DisplayName', 'Sig-spont limit');
        legend('-DynamicLegend')
        axis([OSNRdB(1) OSNRdB(find(OSNRdB ~= 0, 1, 'last')) -8 0])
        xlabel('OSNR (dB)', 'FontSize', 12)
    else
        PrxdBm = Tx.PtxdBm - linkAttdB;
        hline(1) = plot(PrxdBm, log10(ber.count), '-o', 'LineWidth', 2, 'DisplayName', 'Counted');
        hline(2) = plot(PrxdBm, log10(ber.gauss), '-', 'Color', get(hline(1), 'Color'), 'LineWidth', 2, 'DisplayName', 'Gaussian Approximation');
        legend('-DynamicLegend')
        axis([PrxdBm(1) PrxdBm(end) -8 0])
        xlabel('Received power (dB)', 'FontSize', 12)
    end
end
ylabel('log_{10}(BER)', 'FontSize', 12)
set(gca, 'FontSize', 12)