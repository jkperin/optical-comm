function BER = PAM_BER_qsub(fiberLengthKm, wavelengthnm, ModBWGHz, amplified, Ntaps, ENOB, ros)
%% Calculate BER of IM-DD system using M-PAM
% - Equalization is done using a fractionally spaced linear equalizer
% Simulations include modulator, fiber, optical amplifier (optional) characterized 
% only by gain and noise figure, optical bandpass filter, antialiasing 
% filter, sampling, and linear equalization
% Inputs:
% - fiberLength: fiber length in km
% - wavelengthnm: transmission wavelength
% - Modulator: either 'MZM' or 'SiPhotonics'
% - ModBWGHz: modulator bandwidth in GHz. Only used when Modulator == 'SiPhotonics'
% - Amp: whether optical amplifier is included

addpath ../f % general functions
addpath ../apd % for PIN photodetectors

%
ros = eval(ros);
filename = sprintf('results/PAM_BER_L=%skm_lamb=%snm_ModBW=%sGHz_amplified=%s_Ntaps=%s_ENOB=%s_ros=%.2f.mat',...
        fiberLengthKm, wavelengthnm, ModBWGHz, amplified, Ntaps, ENOB, ros);

filename = check_filename(filename)
    
% convert inputs to double (on cluster inputs are passed as strings)
if ~all(isnumeric([fiberLengthKm wavelengthnm ModBWGHz Ntaps, ENOB amplified]))
    fiberLength = 1e3*str2double(fiberLengthKm);
    wavelength = 1e-9*str2double(wavelengthnm);
    ModBW = 1e9*str2double(ModBWGHz);
    Ntaps = round(str2double(Ntaps));
    ENOB = round(str2double(ENOB));
    amplified = logical(str2double(amplified));
end

%% Transmit power swipe
dPower = 0.5; % dBm
if amplified
    Tx.PtxdBm = -23:dPower:-8; % transmitter power range
else
    Tx.PtxdBm = -12:dPower:0; % transmitter power range
end

%% Simulation parameters
sim.Rb = 112e9;    % bit rate in bits/sec
sim.Nsymb = 2^17; % Number of symbols in montecarlo simulation
sim.ros.txDSP = 1; % oversampling ratio transmitter DSP (must be integer). DAC samping rate is sim.ros.txDSP*mpam.Rs
% For DACless simulation must make Tx.dsp.ros = sim.Mct and DAC.resolution = Inf
sim.ros.rxDSP = ros; % oversampling ratio of receiver DSP. If equalization type is fixed time-domain equalizer, then ros = 1
sim.Mct = 2*5;      % Oversampling ratio to simulate continuous time. Must be integer multiple of sim.ros.txDSP and numerator of sim.ros.rxDSP
sim.BERtarget = 1.8e-4; 
sim.Ndiscard = 512; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation
sim.Modulator = 'DML'; % 'MZM' or 'DML'
sim.save = true;
sim.Realizations = 4;
 
%% Simulation control
sim.preAmp = amplified;
sim.preemphasis = false; % preemphasis to compensate for transmitter bandwidth limitation
sim.preemphRange = 25e9; % preemphasis range
sim.mzm_predistortion = 'none'; % predistortion to compensate MZM nonlinearity {'none': no predistortion, 'levels': only PAM levels are predistorted, 'analog': analog waveform is predistorted (DEPRECATED)}
sim.RIN = true; % include RIN noise in montecarlo simulation
sim.phase_noise = true; % whether to simulate laser phase noise
sim.PMD = true; % whether to simulate PMD
sim.quantiz = not(isinf(ENOB)); % whether quantization is included
sim.stopSimWhenBERReaches0 = true; % stop simulation when counted BER reaches 0

% Control what should be plotted
sim.Plots = containers.Map();
sim.Plots('BER') = 0;
sim.Plots('DAC output') = 0;
sim.Plots('Optical eye diagram') = 0;
sim.Plots('Received signal eye diagram') = 0;
sim.Plots('Signal after equalization') = 0;
sim.Plots('Equalizer') = 0;
sim.Plots('Adaptation MSE') = 0;
sim.Plots('Channel frequency response') = 0;
sim.Plots('OSNR') = 0;
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
Tx.DAC.resolution = Inf; % DAC effective resolution in bits
Tx.DAC.rclip = 0;
Tx.DAC.filt = design_filter('bessel', 5, 0.7*mpam.Rs/(sim.fs/2)); % DAC analog frequency response
% juniperDACfit = [-0.0013 0.5846 0];
% Tx.DAC.filt.H = @(f) 1./(10.^(polyval(juniperDACfit, abs(f*sim.fs/1e9))/20)); % DAC analog frequency response

%% Modulator
Tx.rexdB = -15;  % extinction ratio in dB. Defined as Pmin/Pmax
Tx.Mod.type = sim.Modulator; 
Tx.Vgain = 1; % Gain of driving signal
Tx.VbiasAdj = 1; % adjusts modulator bias
Tx.Mod.Vbias = 0.5; % bias voltage normalized by Vpi
Tx.Mod.Vswing = 1;  % normalized voltage swing. 1 means that modulator is driven at full scale
Tx.Mod.BW = ModBW; % DAC frequency response includes DAC + Driver + Modulator
Tx.Mod.filt = design_filter('two-pole', Tx.Mod.BW, sim.fs);
% Tx.Mod.filt = design_filter('bessel', 5, Tx.Mod.BW/(sim.fs/2)); 
Tx.Mod.H = Tx.Mod.filt.H(sim.f/sim.fs);
Tx.Mod.alpha = 0;

%% Laser
% Laser constructor: laser(lambda (nm), PdBm (dBm), RIN (dB/Hz), linewidth (Hz), frequency offset (Hz))
% lambda : wavelength (nm)
% PdBm : output power (dBm)
% RIN : relative intensity noise (dB/Hz)
% linewidth : laser linewidth (Hz)
% freqOffset : frequency offset with respect to wavelength (Hz)
Tx.Laser = laser(wavelength, 0, -150, 0.2e6, 0);
Tx.Laser.H = @(f) Tx.Mod.filt.H(f/sim.fs); % only used if sim.Modulator = 'DML'

%% ============================= Fiber ====================================
% fiber(Length in m, anonymous function for attenuation versus wavelength
% (default: att(lamb) = 0 i.e., no attenuation), anonymous function for 
% dispersion versus wavelength (default: SMF28 with lamb0 = 1310nm, S0 = 0.092 s/(nm^2.km))
Fiber = fiber(fiberLength); 

linkAttdB = Fiber.att(Tx.Laser.wavelength)*Fiber.L/1e3;

%% ========================== Amplifier ===================================
% Constructor: OpticalAmplifier(Operation, param, Fn, Wavelength)
% - Opertation: either 'ConstantOutputPower' or 'ConstantGain'
% - param: GaindB if Operation = 'ConstantGain', or outputPower
% if Operation = 'ConstantOutputPower'
% - Fn:  noise figure in dB
% - Wavelength: operationl wavelength in m
Rx.OptAmp = OpticalAmplifier('ConstantOutputPower', 0, 5, Tx.Laser.wavelength); 
% Rx.OptAmp = OpticalAmplifier('ConstantGain', 20, 5, Tx.Laser.wavelength);

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
Rx.ADC.ENOB = ENOB; % effective number of bits. Quantization is only included if sim.quantiz = true and ENOB ~= Inf
Rx.ADC.rclip = 0;

%% Equalizer
% Terminology: TD = time domain, SR = symbol-rate, LE = linear equalizer
Rx.eq.ros = sim.ros.rxDSP;
Rx.eq.type = 'Adaptive TD-LE';
Rx.eq.Ntaps = Ntaps; % number of taps
Rx.eq.mu = 1e-2; % adaptation ratio
Rx.eq.Ntrain = 1.7e4; % Number of symbols used in training (if Inf all symbols are used)
Rx.eq.Ndiscard = [2e4 128]; % symbols to be discard from begining and end of sequence due to adaptation, filter length, etc
% Note: Ntrain must be at least 4 times larger for ros = 5/4 than for ros
% = 2. Ntrain = 1.5e4 seems sufficient for ros = 5/4

% check if there are enough symbols to perform simulation
Ndiscarded = sum(Rx.eq.Ndiscard) + 2*sim.Ndiscard;
assert(Ndiscarded < sim.Nsymb, 'There arent enough symbols to perform simulation. Nsymb must be increased or Ndiscard must be reduced')
fprintf('%d (2^%.2f) symbols will be used to calculated the BER\n', sim.Nsymb - Ndiscarded, log2(sim.Nsymb - Ndiscarded));

%% Run simulation
%% Calculate BER
Ptx = dBm2Watt(Tx.PtxdBm); % Transmitted power

BER.gauss = zeros(size(Ptx));
BER.count = zeros(size(Ptx));
OSNRdB = zeros(size(Ptx));
for k = 1:length(Ptx)
    Tx.Ptx = Ptx(k);
         
    % Montecarlo simulation
    parfor kk = 1:sim.Realizations
        [BERcount(kk), BERgauss(kk), OSNRdBk(kk)] = ber_pam_montecarlo(mpam, Tx, Fiber, Rx, sim);
    end
    
    BER.count(k) = mean(BERcount);
    BER.gauss(k) = mean(BERgauss);
    OSNRdB(k) = mean(OSNRdBk);
    
    if sim.stopSimWhenBERReaches0 && BER.count(k) == 0
        break;
    end 
end

%% Plots
if sim.shouldPlot('BER') && length(BER.count) > 1
    figure(1), hold on, box on
    if sim.preAmp
        hline(1) = plot(OSNRdB, log10(BER.count), '-o', 'LineWidth', 2, 'DisplayName', 'Counted');
        hline(2) = plot(OSNRdB, log10(BER.gauss), '-', 'Color', get(hline(1), 'Color'), 'LineWidth', 2, 'DisplayName', 'Gaussian Approx.');
        %hline(3) = plot(OSNRdB(OSNRdB ~= 0), log10(pam_ber_from_osnr(mpam.M, OSNRdB(OSNRdB ~= 0), Rx.noiseBW)), '--k', 'LineWidth', 2, 'DisplayName', 'Sig-spont & noise enhancement');
        hline(3) = plot(OSNRdB(OSNRdB ~= 0), log10(pam_ber_from_osnr(mpam.M, OSNRdB(OSNRdB ~= 0), mpam.Rs/2)), ':k', 'LineWidth', 2, 'DisplayName', 'Sig-spont limit');
        legend('-DynamicLegend')
        axis([OSNRdB(1) OSNRdB(find(OSNRdB ~= 0, 1, 'last')) -8 0])
        xlabel('OSNR (dB)', 'FontSize', 12)
    else
        PrxdBm = Tx.PtxdBm - linkAttdB;
        hline(1) = plot(PrxdBm, log10(BER.count), '-o', 'LineWidth', 2, 'DisplayName', 'Counted');
        hline(2) = plot(PrxdBm, log10(BER.gauss), '-', 'Color', get(hline(1), 'Color'), 'LineWidth', 2, 'DisplayName', 'Gaussian Approximation');
        legend('-DynamicLegend')
        axis([PrxdBm(1) PrxdBm(end) -8 0])
        xlabel('Received power (dB)', 'FontSize', 12)
    end
end
ylabel('log_{10}(BER)', 'FontSize', 12)
set(gca, 'FontSize', 12)

%% Save results
if sim.save   
    % delete large variables
    sim = rmfield(sim, 'f');
    sim = rmfield(sim, 't');
    Tx.Mod = rmfield(Tx.Mod, 'H');    
    save(filename)
end