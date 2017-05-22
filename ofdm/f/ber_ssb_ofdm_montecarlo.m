function [ber_count, ber_gauss, SNRndB, OSNRdB] = ber_ssb_ofdm_montecarlo(ofdm, Tx, Fibers, Rx, sim)     
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

%% ================================ DAC ===================================
% Define excursion limits of the DAC
sig = sqrt(0.5*ofdm.var(ofdm.Pn.*abs(Tx.Hdac).^2)); % factor of 0.5 is to account for real and imag parts
if isfield(sim, 'quantiz') && sim.quantiz && not(isinf(Tx.DAC.resolution))
    Tx.DAC.excursion = [-1 1]*Tx.rclip*sig;
end
xt = dac(real(xd), Tx.DAC, sim)... % I
     + 1j*dac(imag(xd), Tx.DAC, sim); % Q

%% ============================= Driver ===================================
xt = xt/(Tx.rclip*sig); % approximately -1 < Re{xt} < 1, -1 < Im{xt} < 1
xt = Tx.Mod.Vswing/2*xt + Tx.Mod.Vbias*(1 + 1j);

%% ============================= Modulator ================================
Tx.Laser.PdBm = Watt2dBm(Tx.Ptx);
Ecw = Tx.Laser.cw(sim);
Etx = mzm(Ecw, xt, Tx.Mod); % transmitted electric field
% Note: Driving signal xt must be normalized by Vpi

% Adjust power to make sure desired power is transmitted
Etx = Etx*sqrt(Tx.Ptx/mean(abs(Etx).^2));

% Carrier-to-signal power ratio
CSPRdB = 10*log10(abs(mean(Etx)).^2/var(Etx))

power_meter(Etx)

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
disp('- IMPORTANT: Simulation including optical amplifier!')
[Erx, OSNRdBtheory, Nase] = Rx.OptAmp.amp(Erx, sim.fs);
linkGain = linkGain*Rx.OptAmp.Gain;

10*log10(mean(abs(Erx(1, :) - mean(Erx(1, :))).^2)/(2*Rx.OptAmp.Ssp*12.5e9))

% Measure OSNR
Osa = OSA(0.3); % optical spectrum analyser with resolution 0.1nm
[~, OSNRdBmeasured] = Osa.estimate_osnr([Erx(1, :) - mean(Erx(1, :)); Erx(2, :) - mean(Erx(2, :))],...
    Tx.Laser.wavelength, sim.f, sim.shouldPlot('OSNR'));
% Note: carrier is removed before measurement

fprintf('OSNR = %.2f dB (theory)\nOSNR = %.2f dB (measured)\n', OSNRdBtheory, OSNRdBmeasured)
OSNRdB = OSNRdBtheory;

%% ========================== Receiver ====================================
%% Direct detection and add thermal noise
% PD.detect(Ein: Input electric field, fs: sampling rate of samples in Ein,
% noise statistics {'gaussian', 'no noise'}, N0: one-sided PSD of thermal
% noise)
fprintf('Photodiode: received power = %.2f dBm\n', power_meter(Erx));
if sim.SSBIcancellation
    % Detection ignoring constant terms, noise-noise beating, and
    % signal-signal beating
    C = mean(Erx(1, :));
    Nase = Nase(1, :);
    S = Erx(1, :) - C - Nase;
    yt = Rx.PD.R*(S.*conj(C + Nase) + C.*conj(S + Nase) + Nase.*conj(S + C)) + sqrt(Rx.N0*sim.fs/2 + Rx.PD.varShot(mean(abs(Erx(1, :)).^2), sim.fs/2))*randn(size(S));
    yt = real(yt);
else
    yt = Rx.PD.detect(Erx, sim.fs, 'gaussian', Rx.N0);
end 

linkGain = linkGain*Rx.PD.Geff;

%% ADC
% ADC performs filtering, quantization, and downsampling
% For an ideal ADC, ADC.ENOB = Inf
if isfield(sim, 'quantiz') && sim.quantiz && not(isinf(Rx.ADC.ENOB))
    sig = std(yt); % = Rx.PD.R*Rx.OptAmp.Gain*sqrt(2*Pc*Ps)
    Rx.ADC.excursion = [-Tx.rclip*sig Tx.rclip*sig];
end

% Rx.ADC.timeRefSignal = abs(xt).^2;
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
