function [ber_count, ber_gauss, OSNRdB, Rx] = ber_pam_montecarlo(mpam, Tx, Fibers, Rx, sim)
%% Calculate BER of M-PAM IM-DD system through montecarlo simulation
% Inputs:
% - mpam: PAM class
% - Tx: struct with transmitter paramters
% - Fibers: array of Fiber classes containing the fibers used in
% transmittion i.e., Fibers = [SMF DCF] or Fibers = [SMF];
% - Rx: struct with receiver parameters
% - sim: struct with simulation parameters

mpamRef = mpam;
mpam = mpam.unbias;

dataTX = randi([0 mpam.M-1], 1, sim.Nsymb); % Random sequence 
dataTX([1 end]) = 0; % set first and last to zero to ensure periodicity

if isfield(sim, 'mzm_predistortion') && strcmpi(sim.mzm_predistortion, 'levels') %% Ordinary PAM with predistorted levels
    % Ajust levels to desired transmitted power and extinction ratio
    mpamPredist = mpam.mzm_predistortion(Tx.Mod.Vswing, Tx.Mod.Vbias, sim.shouldPlot('PAM levels MZM predistortion'));
    xd = mpamPredist.signal(dataTX); % Modulated PAM signal
else
    xd = mpam.signal(dataTX); % Modulated PAM signal
end  

%% ============================ Preemphasis ===============================
if isfield(sim, 'preemphasis') && sim.preemphasis
    femph = abs(freq_time(sim.Nsymb*sim.ros.txDSP, mpam.Rs*sim.ros.txDSP));
    femph(femph >= sim.preemphRange) = 0;
    preemphasis_filter = 10.^(polyval([-0.0013 0.5846 0], femph/1e9)/20);  % Coefficients were measured in the lab  

    xd = real(ifft(fft(xd).*ifftshift(preemphasis_filter)));
end

%% Predistortion to compensate for MZM non-linear response
% This predistorts the analog waveform. If sim.mzm_predistortion ==
% 'levels', then only the levels are predistorted, which is more realistic
if isfield(sim, 'mzm_predistortion') && strcmpi(sim.mzm_predistortion, 'analog')    
    xd = 2/pi*asin(sqrt(abs(xd))).*sign(xd); % apply predistortion
end  

%% ================================ DAC ===================================
% Set DAC time offset in order to remove group delay due to pulse shaping. 
% This way the first sample of xt will be the center of the first pulse. 
% This is only important for plotting.
% Note: Driving signal xd must be normalized by Vpi
Tx.DAC.offset = sim.Mct/mpam.pulse_shape.sps*(length(mpam.pulse_shape.h)-1)/2;
xt = dac(xd, Tx.DAC, sim, sim.shouldPlot('DAC output'));

%% Driver
% Adjust gain to compensate for preemphasis
xt = Tx.Vgain*(xt - mean(xt)) + Tx.VbiasAdj*mean(xt);

%% ============================= Modulator ================================
Tx.Laser.PdBm = Watt2dBm(Tx.Ptx);
Tx.Laser.H = @(f) Tx.Mod.filt.H(f/sim.fs);
if strcmp(Tx.Mod.type, 'MZM')
    Ecw = Tx.Laser.cw(sim);
    Etx = mzm(Ecw, xt, Tx.Mod); % transmitted electric field
elseif strcmp(Tx.Mod.type, 'DML')
    xt = xt - min(xt);
    Etx = Tx.Laser.modulate(xt, sim);
else
    error('ber_pam_montecarlo: Invalid modulator type. Expecting Tx.Mod.type to be either MZM or DML')
end

% Adjust power to ensure that desired power is transmitted
Padj = sqrt(Tx.Ptx/mean(abs(Etx).^2));
if abs(Padj-1) > 0.01
    Etx = Etx*Padj;
    fprintf('Laser output power was adjusted by %.2f to meet target launched power\n', Padj);
end
fprintf('> Launched power: %.2f dBm\n', power_meter(Etx));

% Calculate extinction ratio and intensity levels position
P = abs(Etx).^2;
Psamp = P(1:sim.Mct:end);
Pl = zeros(1, mpam.M);
for k = 1:mpam.M
    Pl(k) = mean(Psamp(dataTX == k-1));
end
Pl = sort(Pl, 'ascend');
rexdB = 10*log10(min(Pl)/max(Pl));
fprintf('Estimated extinction ratio = %.2f dB\n', rexdB)
Pl = Pl - Pl(1);
Pl = (mpam.M-1)*Pl/Pl(end);
fprintf('Normalized optical levels: %s\n', sprintf('%.3f\t', Pl));

%% ========================= Fiber propagation ============================
Erx = Etx;
attdB = 0;
% link_gain = Amp.Gain*Rx.PD.R;
for k = 1:length(Fibers)
    fiberk = Fibers(k); 
    attdB = attdB + fiberk.att(Tx.Laser.wavelength)*fiberk.L/1e3;
    Erx = fiberk.linear_propagation(Erx, sim.f, Tx.Laser.wavelength); % propagation through kth fiber in Fibers
end

%% ========================= Preamplifier =================================
OSNRdB = Inf; % only meaningful when there's a pre-amplifier
if isfield(sim, 'preAmp') && sim.preAmp % only included if sim.preAmp is true
    disp('- IMPORTANT: Simulation includes optical amplifier!')
    [Erx, OSNRdBtheory] = Rx.OptAmp.amp(Erx, sim.fs);
  
    % Measure OSNR
    Osa = OSA(0.3); % optical spectrum analyser with course enough resolution to measure most of signal power
    [~, OSNRdBmeasured] = Osa.estimate_osnr(Erx, Tx.Laser.wavelength, sim.f, sim.shouldPlot('OSNR'));
    % Note: second parameter of Osa.estimate_osnr is OSNR in 0.1 nm
    % resolution
    
    fprintf('OSNR = %.2f dB (theory)\nOSNR = %.2f dB (measured)\n', OSNRdBtheory, OSNRdBmeasured)
    OSNRdB = OSNRdBtheory;
end

%% ========================== Receiver ====================================
%% Direct detection
% Direct detection and add thermal noise
% PD.detect(Ein: Input electric field, fs: sampling rate of samples in Ein,
% noise statistics {'gaussian', 'no noise'}, N0: one-sided PSD of thermal
% noise)
[PrxdBm, Prx] = power_meter(Erx);
yt = Rx.PD.detect(Erx, sim.fs, 'gaussian', Rx.N0);

fprintf('> Received power: %.2f dBm\n', PrxdBm);

% Gain control
yt = yt - mean(yt);
yt = yt*sqrt(mean(abs(mpam.a).^2)/(sqrt(2)*mean(abs(yt).^2)));

%% ADC
% ADC performs filtering, quantization, and downsampling
% For an ideal ADC, ADC.ENOB = Inf
% Rx.ADC.offset = -1; % time offset when sampling      
Hrx = Rx.ADC.filt.H(sim.f/sim.fs);
Rx.ADC.timeRefSignal = xt;
[yk, ~, ytf] = adc(yt, Rx.ADC, sim);

%% =========================== Equalization ===============================
Rx.eq.trainSeq = dataTX; % training sequence
[yd, Rx.eq] = equalize(Rx.eq, yk, [], mpam, sim, sim.shouldPlot('Equalizer'));

% Symbols to be discard in BER calculation
ndiscard = [1:Rx.eq.Ndiscard(1)+sim.Ndiscard (sim.Nsymb-Rx.eq.Ndiscard(2)-sim.Ndiscard):sim.Nsymb];
ydfull = yd;
yd(ndiscard) = []; 
dataTX(ndiscard) = [];

%% =========================== Detection ==============================
% dataRX = mpam.demod(yd);
[dataRX, mpam] = mpam.demod_sweeping_thresholds(yd, dataTX);
% Note: when performing demodulaton sweepingn thresholds, the decision
% thresholds will be picked to minimize BER

%% Counted BER
[~, ber_count] = biterr(dataRX, dataTX);

%% ====================== Gaussian approximation ===========================
% Note 1: Gaussian approximation is used despite the fact that the dominant
% noise in pre-amplified systems is signal-spontaneous beat noise, which is
% not Gaussian distributed.
% Note 2: This calculation includes noise enhancement due to equalization.

% Noise bandwidth
Hrxeq = abs(Hrx.*Rx.eq.Hff(sim.f/(Rx.eq.ros*mpam.Rs))).^2;
Hrxeq = Hrxeq/interp1(sim.f, Hrxeq, 0); % normalize to have unit gain at DC
noiseBW = 0.5*trapz(sim.f, Hrxeq);
Rx.noiseBW = noiseBW;

if isfield(sim, 'preAmp') && sim.preAmp % amplified system: signal-spontaneous beat noise dominant
    BWopt = sim.fs; % optical filter (OF) noise bandwidth = sampling rate, since OF is not included 
    % Not divided by 2 because optical filter is a bandpass filter

    % Noise std for intensity level Plevel
    Npol = 2; % number of polarizations. Npol = 1, if polarizer is present, Npol = 2 otherwise.
    
    noiseSTD = @(Plevel) sqrt(noiseBW*Rx.N0 + Rx.PD.varShot(Plevel, noiseBW)... % thermal + shot
                    + Rx.PD.R^2*(Rx.OptAmp.varNoiseDD(Plevel/(Rx.OptAmp.Gain), noiseBW, BWopt, Npol))); % sig-spont + spont-spont    
    % Note 1: Plevel is divided by amplifier gain to obtain power at the amplifier input
    % Note 2: if the amplifier is operating in the constant output power 
    % mode, the appropriate amplifier gain was already set in montecarlo
    % simulation. The BER estimation function will assume the same amplifier gain.

else % unamplified system: thermal-noise dominant
    noiseSTD = Rx.PD.stdNoise(Hrx, Rx.eq.Hff(sim.f/(Rx.eq.ros*mpam.Rs)), Rx.N0, Tx.Laser.RIN, sim);
end

% AWGN approximation
mpamRef = mpamRef.adjust_levels(Prx, -Inf); % may replace -Inf for rexdB
ber_gauss = mpamRef.berAWGN(noiseSTD);

%% Plots
if sim.shouldPlot('Optical eye diagram')  
    Nstart = sim.Ndiscard*sim.Mct + 1;
    figure(103), clf, hold on, box on
    eyediagram(abs(Etx(Nstart:end)).^2, 2*sim.Mct)
    a = axis;
    plot(a(1:2), [1; 1]*Pl, '-k', 'Linewidth', 2)
    title('Transmitted optical signal eye diagram')
    drawnow
end

if sim.shouldPlot('Received signal eye diagram')
    mpam = mpam.norm_levels();
    Nstart = sim.Ndiscard*sim.Mct + 1;
    figure(104), clf, box on, hold on
    eyediagram(ytf(Nstart:end), 2*sim.Mct)
    title('Received signal eye diagram')
    a = axis;
    h1 = plot(a(1:2), (mpam.a*[1 1]).', '-k');
    h2 = plot(a(1:2), (mpam.b*[1 1]).', '--k');
    h3 = plot((sim.Mct+1)*[1 1], a(3:4), 'k');
    legend([h1(1) h2(1) h3], {'Levels', 'Decision thresholds', 'Sampling point'})
    drawnow
end

if sim.shouldPlot('Signal after equalization')
    mpam = mpam.norm_levels();
    figure(105), clf, box on, hold on
    h1 = plot(ydfull, 'o');
    a = axis;
    h2= plot(a(1:2), (mpam.a*[1 1]).', '-k');
    h3 = plot(a(1:2), (mpam.b*[1 1]).', '--k');
    h4 = plot((Rx.eq.Ndiscard(1)+sim.Ndiscard)*[1 1], a(3:4), ':k');
    plot(((sim.Nsymb-Rx.eq.Ndiscard(2)-sim.Ndiscard))*[1 1], a(3:4), ':k');
    legend([h1 h2(1) h3(1) h4], {'Equalized samples', 'PAM levels',...
        'Decision thresholds', 'BER measurement window'})
    title('Signal after equalization')
%     axis([1 sim.Nsymb -0.2 1.2])
    drawnow
end

if sim.shouldPlot('Heuristic noise pdf')
    figure(106)
    [nn, xx] = hist(yd(dataTX == 2), 50);
    nn = nn/trapz(xx, nn);
    bar(xx, nn)
    drawnow
end

if sim.shouldPlot('Decision errors')
    figure(107)
    stem(dataRX ~= dataTX)
    ylabel('Errors')
    xlabel('Symbol')
    title('Decision errors')
    drawnow
end