function [PrxdBm_BERtarget_opt, PrxdBm_BERtarget_pin, PrxdBm_BERtarget]...
        = sensitivity_vs_L_qsub(M, ka, level_spacing, BW0GHz, GainBWGHz, modBWGHz, lamb, Lkm)
%% Function to calculate receiver sensitivity for a specified fiber length 
% On cluster(corn) inputs are passed as strings
% Inputs:
% - M = PAM order
% - ka = impact ionization factor
% - level_spacing = {eqully-spaced, optimized}
% - BW0GHz = low-gain APD bandwidth in GHz
% - GainBWGHz = Gain-bandwidth product in GHz
% - modBW = modulator BW in GHz
% Note: BW0GHz = Inf and GainBWGHz = Inf and modBW = Inf, this corresponds
% to AWGN simulation
% Outputs:
% - PrxdBm_BERtarget_opt: received power necessary to achieve target BER
% for APD with optimal gain
% - PrxdBm_BERtarget_pin: received power necessary to achieve target BER
% for PIN receiver with infinite bandwidth 
% - PrxdBm_BERtarget: received power necessary to achieve target BER
% with APD gain set for several values near the optimal gain (used for
% verification)

addpath ../mpam
addpath ../f
addpath f

% Convert inputs to double
if ~isnumeric([M, ka, BW0GHz, GainBWGHz, modBWGHz])
    M = round(str2double(M));
    ka = str2double(ka);
    BW0GHz = str2double(BW0GHz);
    GainBWGHz = str2double(GainBWGHz);
    modBWGHz = str2double(modBWGHz);
    lamb = str2double(lamb);
    Lkm = str2double(Lkm);
end

filename = sprintf('results/sensitivity_vs_L/sensitivity_vs_L_%d-PAM_%s_ka=%d_BW0=%d_GainBW=%d_modBW=%d_lamb=%dnm_L=%dkm',...
    M, level_spacing, round(100*ka), BW0GHz, GainBWGHz, modBWGHz, lamb, Lkm)

%% Simulation parameters
sim.Nsymb = 2^18; % Number of symbols in montecarlo simulation
sim.ros.txDSP = 1; % oversampling ratio of transmitter DSP
sim.ros.rxDSP = 1; % oversampling ratio of receivre DSP
sim.Mct = 8*sim.ros.rxDSP;    % Oversampling ratio to simulate continuous time (must be odd so that sampling is done right, and FIR filters have interger grpdelay)  
sim.L = 4;        % de Bruijin sub-sequence length (ISI symbol length)
sim.BERtarget = 1.8e-4; 
sim.Ndiscard = 1024;  % number of 0 symbols to be inserted at the begining and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

%
sim.WhiteningFilter = true;
sim.save = true;
sim.RIN = true; % include RIN noise. Only included in montecarlo simulation
sim.quantiz = false; % whether quantization is included
sim.terminateWhenBERReaches0 = true; % whether simulation is terminated when counted BER reaches 0

% What to plot
sim.Plots = containers.Map();
% sim.Plots('Gain swipe') = 1;
% sim.Plots('BER') = 1;
% sim.Plots('Adaptation MSE') = 0;
% sim.Plots('Frequency Response') = 0;
% sim.Plots('Equalizer') = 1;
% sim.Plots('DAC output') = 1;
% sim.Plots('Optical eye diagram') = 1;
% sim.Plots('Received signal eye diagram') = 1;
% sim.Plots('Signal after equalization') = 1;
sim.shouldPlot = @(x) sim.Plots.isKey(x) && sim.Plots(x);

%% M-PAM
pulse_shape = select_pulse_shape('rect', sim.ros.txDSP);
mpam = PAM(M, 107e9, level_spacing, pulse_shape);

%% Time and frequency
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'
[sim.f, sim.t] = freq_time(sim.N, sim.fs);

%% Laser
% Laser constructor: laser(lambda (nm), PdBm (dBm), RIN (dB/Hz), linewidth (Hz), frequency offset (Hz))
% lambda : wavelength (nm)
% PdBm : output power (dBm)
% RIN : relative intensity noise (dB/Hz)
% linewidth : laser linewidth (Hz)
% freqOffset : frequency offset with respect to wavelength (Hz)
Tx.Laser = laser(lamb*1e-9, 0, -150, 0.2e6, 0);

%% Transmitter
Tx.PtxdBm = -25:-3; % transmitted power swipe
% Tx.PtxdBm = -15; % transmitted power swipe

%% DAC
Tx.DAC.fs = sim.ros.txDSP*mpam.Rs; % DAC sampling rate
Tx.DAC.ros = sim.ros.txDSP; % oversampling ratio of transmitter DSP
Tx.DAC.resolution = Inf; % DAC effective resolution in bits
Tx.DAC.filt = design_filter('bessel', 5, 100e9/(sim.fs/2)); % DAC analog frequency response

%% Modulator
Tx.Mod.alpha = 2; % chirp parameter
Tx.Mod.rexdB = -10;  % extinction ratio in dB. Defined as Pmin/Pmax
% Modulator frequency response
Tx.Mod.BW = modBWGHz*1e9;
Tx.Mod.fc = Tx.Mod.BW/sqrt(sqrt(2)-1); % converts to relaxation frequency
Tx.Mod.h = @(t) (2*pi*Tx.Mod.fc)^2*t(t >= 0).*exp(-2*pi*Tx.Mod.fc*t(t >= 0));
Tx.Mod.grpdelay = 2/(2*pi*Tx.Mod.fc);  % group delay of second-order filter in seconds
Tx.Mod.H = @(f) exp(1j*2*pi*f.*Tx.Mod.grpdelay)./(1 + 2*1j*f/Tx.Mod.fc - (f/Tx.Mod.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)

%% Fiber
Fiber = fiber(Lkm*1e3);

%% Receiver
Rx.N0 = (30e-12).^2; % thermal noise psd
Rx.Id = 10e-9;
Rx.R = 1;
Rx.filtering = 'matched'; % Electric Lowpass Filter prior to sampling and equalization: either 'antialisaing' or 'matched'

%% ADC
Rx.ADC.ros = sim.ros.rxDSP;
Rx.ADC.fs = Rx.ADC.ros*mpam.Rs*(mpam.pulse_shape.rolloff + 1)/2;
Rx.ADC.filt = design_filter('butter', 5, 0.5*Rx.ADC.fs/(sim.fs/2)); % Antialiasing filter
Rx.ADC.ENOB = Inf; % effective number of bits. Quantization is only included if sim.quantiz = true and ENOB ~= Inf
Rx.ADC.rclip = 0;

%% Equalization
if isinf(BW0GHz) && isinf(modBWGHz) % awgn simulation
    disp('> Equalizer = None')
    Rx.eq.type = 'None';
    Rx.eq.ros = 1;
    Rx.filtering = 'matched'; % Electric Lowpass Filter prior to sampling and equalization: either 'antialisaing' or 'matched'
    Rx.eq.Ndiscard = 1024*[1 1]; % symbols to be discard from begining and end of sequence due to adaptation, filter length, etc
else
    disp('> Equalizer = fixed time-domain symbol-rate linear equalizer')
    Rx.eq.type = 'Fixed TD-SR-LE';
    Rx.eq.ros = sim.ros.rxDSP;
    Rx.eq.Ntaps = 31;
    Rx.eq.Ndiscard = 1024*[1 1]; % symbols to be discard from begining and end of sequence due to adaptation, filter length, etc
    
    % Adaptive
%     disp('> Equalizer = adaptive time-domain symbol-rate linear equalizer')
%     Rx.eq.type = 'Adaptive TD-LE';
%     Rx.eq.ros = sim.ros.rxDSP;
%     Rx.eq.Ntaps = 31;
%     Rx.eq.mu = 1e-2;
%     Rx.eq.Ntrain = 0.7e4; % Number of symbols used in training (if Inf all symbols are used)
%     Rx.eq.Ndiscard = [1e4 512]; % symbols to be discard from begining and end of sequence due to adaptation, filter length, etc
end

%% APD 
% apd(GaindB, ka, BW, R, Id) 
Apd = apd(10, ka, 1e9*[BW0GHz GainBWGHz], Rx.R, Rx.Id);

% Received power 
[~, link_att_dB] = Fiber.link_attenuation(Tx.Laser.wavelength);
PrxdBm = Tx.PtxdBm - link_att_dB;

%% Calculate BER for PIN with same properties, but inifinite bandwidth
PIN = pin(Rx.R, Rx.Id);

ber_pin = apd_ber(mpam, Tx, Fiber, PIN, Rx, sim);

PrxdBm_BERtarget_pin = interp1(log10(ber_pin.enum), PrxdBm, log10(sim.BERtarget));

%% APD optimal gainn
% Find gain for optimal margin
sim.OptimizeGain = true;

% Optimal gain
[ber_apd_opt, ~, Apd] = apd_ber(mpam, Tx, Fiber, Apd, Rx, sim);
Gopt = Apd.Gain;

PrxdBm_BERtarget_opt = interp1(log10(ber_apd_opt.enum), PrxdBm, log10(sim.BERtarget));

%% Swipe gain
sim.OptimizeGain = false;
Gains = Gopt + [-3:-1 1:3];
PrxdBm_BERtarget = zeros(size(Gains));

if sim.shouldPlot('Gain swipe')
    figure(10), clf, hold on, box on
    leg = {};
    hline = [];
    hline(1) = plot(PrxdBm, log10(ber_pin.enum), '-b');
    plot(PrxdBm, log10(ber_pin.awgn), '--', 'Color', get(hline(1), 'Color'))
    plot(PrxdBm, log10(ber_pin.count), '-o', 'Color', get(hline(1), 'Color'))
    leg = [leg 'PIN'];
    
    hline(2) = plot(PrxdBm, log10(ber_apd_opt.enum), '-k');
    plot(PrxdBm, log10(ber_apd_opt.awgn), '--', 'Color', get(hline(2), 'Color'))
    plot(PrxdBm, log10(ber_apd_opt.count), '-o', 'Color', get(hline(2), 'Color'))
    leg = [leg sprintf('Optimal Gain = %.2f', Gopt)];
    
    xlabel('Received Power (dBm)')
    ylabel('log(BER)') 
    axis([PrxdBm(1) PrxdBm(end) -8 0])
    set(gca, 'xtick', PrxdBm)
    title(sprintf('%d-PAM, %s, ka = %.2f, BW0 = %.2f, GainBW = %.2f', mpam.M, mpam.level_spacing, ka, BW0GHz, GainBWGHz))
    drawnow
end

fprintf('- %d-PAM, %s, ka = %.2f, BW0 = %.2f, GainBW = %.2f, L = %.2f\n',...
    mpam.M, mpam.level_spacing, ka, Apd.BW0/1e9, Apd.GainBW/1e9, Lkm)
for k= 1:length(Gains)
    fprintf('-- Gain = %.2f\n', Gains(k))

    Apd.Gain = Gains(k);

    % BER
    ber_apd = apd_ber(mpam, Tx, Fiber, Apd, Rx, sim);
    BER(k) = ber_apd;

    % Calculate power at the target BER
    PrxdBm_BERtarget(k) = interp1(log10(ber_apd.enum), PrxdBm, log10(sim.BERtarget));

    % Figures
    if sim.shouldPlot('Gain swipe')
        hline(k+2) = plot(PrxdBm, log10(ber_apd.enum), '-');
        plot(PrxdBm, log10(ber_apd.awgn), '--', 'Color', get(hline(k+2), 'Color'))
        plot(PrxdBm, log10(ber_apd.count), '-o', 'Color', get(hline(k+2), 'Color'))
        leg = [leg sprintf('Gain = %.2f', Apd.Gain)];
        drawnow
    end
end   
if sim.shouldPlot('Gain swipe')
    legend(hline, leg)  
end

% Save results
if sim.save
    sim = rmfield(sim, 't');
    sim = rmfield(sim, 'f');
    save(filename)
end
