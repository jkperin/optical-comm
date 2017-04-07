function [ber_count, ber_gauss, SNRndB, OSNRdB] = ber_ofdm_montecarlo(ofdm, Tx, Fibers, Rx, sim)     
%% Calculate BER of pre-amplified IM-DD system through montecarlo simulation
% Inputs:
% - ofdm: OFDM class
% - Tx: struct with transmitter paramters
% - Fibers: array of Fiber classes containing the fibers used in
% transmittion i.e., Fibers = [SMF DCF] or Fibers = [SMF];
% - Rx: struct with receiver parameters
% - sim: struct with simulation parameters  

%% Generate OFDM signal
[xd, Rx.AdEq.trainSeq] = ofdm.signal(sim.Nsymb); 

if isfield(sim, 'clipping_compensation') && sim.clipping_compensation
    rclip = Tx.rclip;
    K = 1 - 2*qfunc(rclip);
    sig = std(xd); %sqrt(ofdm.var(ofdm.Pn.*abs(Tx.Hdac).^2));
    xc = xd;
    xc(xd > rclip*sig) = rclip*sig;
    xc(xd < -rclip*sig) = -rclip*sig;
    d = xc - K*xd;
    xd = xd - d/K;
%     xd = xd/K;
end

%% ================================ DAC ===================================
% Define excursion limits of the DAC
if isfield(sim, 'quantiz') && sim.quantiz && not(isinf(Tx.DAC.resolution))
    sig = sqrt(ofdm.var(ofdm.Pn.*abs(Tx.Hdac).^2));
    if strcmpi(ofdm.prefix, 'ACO') % ACO-OFDM
        Tx.DAC.excursion = [0 Tx.rclip*sig];
    else % DC-OFDM
        Tx.DAC.excursion = [-1 1]*Tx.rclip*sig;
    end
end
xt = dac(xd, Tx.DAC, sim); % digital-to-analog conversion 

%% ============================= Driver ===================================
[Padj, dc] = ofdm.adjust_power_allocation(ofdm.Pn.*abs(Tx.Hdac).^2, Tx.Ptx, Tx.rclip, Tx.rexdB); % Factor to refer subcarriers powers to the transmitter
xt = xt*sqrt(Padj) + dc;
ofdm.Pn = ofdm.Pn*Padj;
xt(xt < 0) = 0; % Clip negative excursion

% % Discard first and last symbols
% xt(1:sim.Mct*sim.Ndiscard) = 0; % zero sim.Ndiscard first symbols
% xt(end-sim.Mct*sim.Ndiscard+1:end) = 0; % zero sim.Ndiscard last symbols

%% ============================= Modulator ================================
Tx.Laser.PdBm = Watt2dBm(Tx.Ptx);
Tx.Laser.H = @(f) Tx.Mod.filt.H(f/sim.fs);
if strcmp(Tx.Mod.type, 'MZM')
    xt = Tx.Mod.Vswing*xt/(Vbias + Tx.rclip*sqrt(ofdm.var(Pnnorm))); 
    % Note: Driving signal xd must be normalized by Vpi
    
    Ecw = Tx.Laser.cw(sim);
    Etx = mzm(Ecw, xt, Tx.Mod); % transmitted electric field
elseif strcmp(Tx.Mod.type, 'DML')
    Etx = Tx.Laser.modulate(xt, sim);
    % Note: xt does not need to be normalized here. Normalization is
    % performed later to ensure that desired transmit power is reached
else
    error('ber_pam_montecarlo: Invalid modulator type. Expecting Tx.Mod.type to be either MZM or DML')
end

% Adjust power to make sure desired power is transmitted
Etx = Etx*sqrt(Tx.Ptx/mean(abs(Etx).^2));

%% ========================= Fiber propagation ============================
Erx = Etx;
attdB =0;
linkGain = 1;
for k = 1:length(Fibers)
    fiberk = Fibers(k); 
    attdB = attdB + fiberk.att(Tx.Laser.wavelength)*fiberk.L/1e3;
    linkGain = linkGain*fiberk.link_attenuation(Tx.Laser.wavelength);
    Erx = fiberk.linear_propagation(Erx, sim.f, Tx.Laser.wavelength); % propagation through kth fiber in Fibers
end

%% ========================= Preamplifier =================================
OSNRdB = Inf; % only meaningful when there's a pre-amplifier
if isfield(sim, 'preAmp') && sim.preAmp % only included if sim.preAmp is true
    disp('- IMPORTANT: Simulation including optical amplifier!')
    [Erx, OSNRdBtheory] = Rx.OptAmp.amp(Erx, sim.fs);
    linkGain = linkGain*Rx.OptAmp.Gain;
   
    % Measure OSNR
    Osa = OSA(0.1); % optical spectrum analyser with resolution 0.1nm
    OSNRdBmeasured = Osa.estimate_osnr(Erx, Tx.Laser.wavelength, sim.f, sim.shouldPlot('OSNR'));
    
    fprintf('OSNR = %.2f dB (theory)\nOSNR = %.2f dB (measured)\n', OSNRdBtheory, OSNRdBmeasured)
    OSNRdB = OSNRdBtheory;
end

%% ========================== Receiver ====================================
%% Direct detection and add thermal noise
% PD.detect(Ein: Input electric field, fs: sampling rate of samples in Ein,
% noise statistics {'gaussian', 'no noise'}, N0: one-sided PSD of thermal
% noise)
fprintf('Photodiode: received power = %.2f dBm\n', power_meter(Erx));
yt = Rx.PD.detect(Erx, sim.fs, 'gaussian', Rx.N0);
linkGain = linkGain*Rx.PD.Geff;

%% ADC
% ADC performs filtering, quantization, and downsampling
% For an ideal ADC, ADC.ENOB = Inf
if isfield(sim, 'quantiz') && sim.quantiz && not(isinf(Rx.ADC.ENOB))
    if strcmpi(ofdm.prefix, 'ACO') % ACO-OFDM
        sig = sqrt(ofdm.var(ofdm.Pn*linkGain^2.*abs(sim.Hch).^2));
        Rx.ADC.excursion = [0, Tx.rclip*sig];
    else % DC-OFDM
        sig = std(yt);
        Rx.ADC.excursion = [0, 2*Tx.rclip*sig];
    end
end

Rx.ADC.timeRefSignal = xt;
yk = adc(yt, Rx.ADC, sim);

%% OFDM detection
[Xn, AGCn, W] = ofdm.detect(yk, Rx.AdEq, sim.shouldPlot('Adaptation MSE') || sim.shouldPlot('Constellations') || sim.shouldPlot('Equalizer'));

%% Calculate BER
[ber_count, ~] = ofdm.count_ber([Rx.AdEq.Ntrain+sim.Ndiscard sim.Ndiscard], sim.shouldPlot('Decision errors'));

%% Gaussian approximation
valid_window = Rx.AdEq.Ntrain+sim.Ndiscard:sim.Nsymb-sim.Ndiscard; % valid window for measurements (discard adapation transients)
Xnvalid = Xn(:, valid_window);
Pnrx = mean(abs(Xnvalid).^2, 2)./abs(AGCn.*W).^2; % received power at each subcarrier
ofdm.Pn = ofdm.Pn*linkGain^2; % refer subcarrier powers to the receiver
% Note: the linkGain is an optical power gain. Following from the Gausssian
% approximation of the OFDM signal, the subcarriers powers are scaled by
% linkGain^2.
noise = Xnvalid - Rx.AdEq.trainSeq(:, valid_window);
Pnoise = mean(abs(noise).^2, 2);
Pnoise = Pnoise./abs(AGCn.*W).^2;

[ber_gauss, SNRndB] =  ofdm.estimate_ber(Pnrx.', Pnoise, sim.Hch, sim.varNoise, sim.shouldPlot('Estimated SNR'));
% Note: estimated noise variance in each subcarrier sim.varNoise does not
% incluse noise enhancement
