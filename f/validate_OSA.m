%% Validate OSA OSNR estimate
clear, clc, close all

%% ======================== Simulation parameters =========================
sim.Nsymb = 2^12; % Number of symbols in montecarlo simulation
sim.Mct = 15;    % Oversampling ratio to simulate continuous time 
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation
sim.Rb = 2*112e9; % Bit rate
sim.Npol = 2;                                                              % number of polarizations
sim.pulse_shape = select_pulse_shape('rect', sim.Mct);                     % pulse shape
ModFormat = QAM(4, sim.Rb/sim.Npol, sim.pulse_shape);                  % M-QAM modulation format
% ModFormat = DPSK(4, sim.Rb/sim.Npol, sim.pulse_shape);                     % M-DPSK modulation format

PdBm = -20;

sim.fs = ModFormat.Rs*sim.Mct;
[sim.f, sim.t] = freq_time(sim.N, sim.fs);

filt = design_filter('bessel', 5, 0.7*ModFormat.Rs/(sim.fs/2));
Mod.H = filt.H(sim.f/sim.fs);

%% Generate transmitted symbols
dataTX = randi([0 ModFormat.M-1], [2, sim.Nsymb]); % symbol stream for each polarization
Nzero = 10; % zero Nzero first and last symbols to make sequence periodic
dataTX(:, 1:Nzero) = 0;
dataTX(:, end-Nzero+1:end) = 0;

[Vin, symbolsTX] = ModFormat.signal(dataTX); % X & Y pol

% Filter drive waveforms for modulators txfilt.
% group delay of Tx.filt.H has already been removed
Htx = ifftshift(filt.H(sim.f/sim.fs).*exp(1j*2*pi*sim.f/sim.fs*ModFormat.pulse_shape_grpdelay)); % transmitter filter and remove group delay due to pulse shaping in ModFormat
Vout(1, :) = real(ifft(fft(real(Vin(1, :))).*Htx)) + 1j*real(ifft(fft(imag(Vin(1, :))).*Htx)); 
Vout(2, :)= real(ifft(fft(real(Vin(2, :))).*Htx)) + 1j*real(ifft(fft(imag(Vin(2, :))).*Htx));

Laser = laser(1550e-9, PdBm, -150, 200);

Ein = Laser.cw(sim); % Generates electric field with intensity and phase noise       
Ein = mzm(Ein, Vout, Mod); % modulate optical signal using eletro-optical modulator (EOM)

% Ensure that transmitted power is a desired level
% Ein(2, :) = 0;
Ein = Ein*sqrt(dBm2Watt(PdBm)/sum(mean(abs(Ein).^2, 2)));

fprintf('Optical power before amplifier = %.2f dBm\n', power_meter(Ein))

%% ======================== Optical Amplifier =============================
% Constructor: OpticalAmplifier(Operation, param, Fn, Wavelength)
% - Opertation: either 'ConstantOutputPower' or 'ConstantGain'
% - param: GaindB if Operation = 'ConstantGain', or outputPower
% if Operation = 'ConstantOutputPower'
% - Fn:  noise figure in dB
% - Wavelength: operationl wavelength in m
OptAmp = OpticalAmplifier('ConstantOutputPower', 0, 5, Laser.wavelength);
% Rx.OptAmp = OpticalAmplifier('ConstantGain', 20, 5, Tx.Laser.wavelength);
% Note: the amplifier here operates in the constant output power mode,
% where the output power after amplification is set to Rx.AmpOutPowerdBm

[Eout, OSNRdBtheory, W] = OptAmp.amp(Ein, sim.fs);

[~, Pout] = power_meter(Eout);
[~, Pnoise] = power_meter(W);
OSNRdBcheck = 10*log10((Pout/Pnoise-1)*sim.fs/(2*12.5e9))


% Measure OSNR
Osa = OSA(0.5); % optical spectrum analyser with resolution 0.1nm
[OSNRdBmeasured, OSNRdBmeasured01nm] = Osa.estimate_osnr(Eout, Laser.wavelength, sim.f, true);

fprintf('OSNR = %.2f dB (theory)\nOSNR = %.2f dB (measured)\n', OSNRdBtheory, OSNRdBmeasured)


