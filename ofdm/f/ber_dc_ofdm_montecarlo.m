function [ber_count, ber_gauss, SNRndB, OSNRdB] = ber_dc_ofdm_montecarlo(ofdm, Tx, Fibers, Rx, sim)     
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
xd = xd/sqrt(ofdm.var); % normalize signal, so that var(xd) = 1
Pnnorm =ofdm.Pn/ofdm.var;
Vbias = ofdm.dc_bias(Pnnorm, Tx.rclip, Tx.Mod.Hdac, Tx.rexdB);
xd = xd + Vbias;

%% ================================ DAC ===================================
Tx.DAC.excursion = [0 Vbias+Tx.rclip*sqrt(ofdm.var(Pnnorm))];
xt = dac(xd, Tx.DAC, sim);

% % Discard first and last symbols
% xt(1:sim.Mct*sim.Ndiscard) = 0; % zero sim.Ndiscard first symbols
% xt(end-sim.Mct*sim.Ndiscard+1:end) = 0; % zero sim.Ndiscard last symbols

%% ============================= Modulator ================================
%% Driver
xt(xt < 0) = 0; % DC bias was already ajusted. So just clip at zero

Tx.Laser.PdBm = Watt2dBm(Tx.Ptx);
Tx.Laser.H = @(f) Tx.Mod.filt.H(f/sim.fs);
if strcmp(Tx.Mod.type, 'MZM')
    xt = xt/(Vbias*Tx.rclip); % normalize so that excursion is [0, 1] (full swing)
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
% link_gain = Amp.Gain*Rx.PD.R;
for k = 1:length(Fibers)
    fiberk = Fibers(k); 
    attdB = attdB + fiberk.att(Tx.Laser.wavelength)*fiberk.L/1e3;
    Erx = fiberk.linear_propagation(Erx, sim.f, Tx.Laser.wavelength); % propagation through kth fiber in Fibers
end

%% ========================= Preamplifier =================================
OSNRdB = Inf; % only meaningful when there's a pre-amplifier
if isfield(sim, 'preAmp') && sim.preAmp % only included if sim.preAmp is true
    disp('- IMPORTANT: Simulation including optical amplifier!')
    [Erx, OSNRdBtheory] = Rx.OptAmp.amp(Erx, sim.fs);
   
    % Adjust power to pre-defined value
    Att = dBm2Watt(Rx.OptAmpOutPowerdBm)/dBm2Watt(power_meter(Erx));
    Erx = Erx*sqrt(Att);  % keep constant received power

    % Measure OSNR
    Osa = OSA(0.1); % optical spectrum analyser with resolution 0.1nm
    OSNRdBmeasured = Osa.estimate_osnr(Erx, Tx.Laser.wavelength, sim.f, sim.shouldPlot('OSNR'));
    
    fprintf('OSNR = %.2f dB (theory)\nOSNR = %.2f dB (measured)\n', OSNRdBtheory, OSNRdBmeasured)
    OSNRdB = OSNRdBtheory;
    % check
%         OSNRdBtheorycheck = 10*log10(dBm2Watt(Tx.PlaunchdBm(k))/(2*Rx.OptAmp.nsp*Rx.OptAmp.h*Rx.OptAmp.c/Tx.Laser.lambda*12.5e9))
%         SNRdBtheorycheck = 10*log10(dBm2Watt(Tx.PlaunchdBm(k))/(Rx.OptAmp.nsp*Rx.OptAmp.h*Rx.OptAmp.c/Tx.Laser.lambda*sim.Rs))
end

%% ========================== Receiver ====================================
%% Direct detection and add thermal noise
% PD.detect(Ein: Input electric field, fs: sampling rate of samples in Ein,
% noise statistics {'gaussian', 'no noise'}, N0: one-sided PSD of thermal
% noise)
[PrxdBm, Prx] = power_meter(Erx);
yt = Rx.PD.detect(Erx, sim.fs, 'gaussian', Rx.N0);
yref = yt;

fprintf('> Received power: %.2f dBm\n', PrxdBm);

%% Automatic gain control
yt = yt - mean(yt);
yt = yt/mean(abs(yt).^2);

%% ADC
% ADC performs filtering, quantization, and downsampling
% For an ideal ADC, ADC.ENOB = Inf
yk = adc(yt, Rx.ADC, sim);

%% OFDM detection
Xn = ofdm.detect(yk, Rx.AdEq, sim.shouldPlot('Adaptation MSE') || sim.shouldPlot('Constellations') || sim.shouldPlot('Equalizer'));

%% Calculate BER
[ber_count, ~] = ofdm.countBER([Rx.AdEq.Ntrain+sim.Ndiscard sim.Ndiscard]);

%% AWGN approximation
[ber_gauss, SNRndB] =  ofdm.estimate_ber(yref, sim.Hch, sim.varNoise, sim.shouldPlot('Estimated SNR'));
