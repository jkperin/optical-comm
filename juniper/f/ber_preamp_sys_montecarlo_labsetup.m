function [ber_count, ber_awgn, OSNRdB] = ber_preamp_sys_montecarlo_labsetup(mpam, Tx, Fibers, Amp, Rx, sim)
%% Calculate BER of pre-amplified IM-DD system through montecarlo simulation
% Inputs:
% - mpam: PAM class
% - Tx: struct with transmitter paramters
% - Fibers: array of Fiber classes containing the fibers used in
% transmittion i.e., Fibers = [SMF DCF] or Fibers = [SMF];
% - Amp: pre-amplifier using SOA class
% - Rx: struct with receiver parameters
% - sim: struct with simulation parameters

mpamRef = mpam;
mpam = mpam.unbias;

dataTX = randi([0 mpam.M-1], 1, sim.Nsymb); % Random sequence 
Nzero = 10;
dataTX([1:Nzero end-Nzero+1:end]) = 0; % set first and last Nzero symbols to 0 to make sequence periodic

%% ======================= Duobinary enconding  ===========================
if isfield(sim, 'duobinary') && sim.duobinary
    mpamdb = mpam.set_levels(0:mpam.M-1, 0.5 + (0:mpam.M-2));    
    
    xd = mpamdb.signal(dataTX); % Modulated PAM signal
    
    xd = duobinary_encoding(xd);
    
    xd_enc = xd;
else %% Ordinary PAM
    % Ajust levels to desired transmitted power and extinction ratio
    if isfield(sim, 'mzm_predistortion') && strcmpi(sim.mzm_predistortion, 'levels')   
        mpamPredist = mpam.mzm_predistortion(Tx.Mod.Vswing, Tx.Mod.Vbias, sim.shouldPlot('PAM levels MZM predistortion'));
        mpamPredist = mpamPredist.unbias;
        xd = mpamPredist.signal(dataTX); % Modulated PAM signal
    else
        xd = mpam.signal(dataTX); % Modulated PAM signal
    end 
end  

%% ============================ Preemphasis ===============================
if sim.preemphasis
    femph = abs(freq_time(sim.Nsymb*sim.ros.txDSP, mpam.Rs*sim.ros.txDSP));
    femph(femph >= sim.preemphRange) = 0;
    emphasis_filter = 10.^(polyval([-0.0013 0.5846 1.5859], femph/1e9)/20);    

    xd = real(ifft(fft(xd).*ifftshift(emphasis_filter)));
end

%% ============================== Driver ==================================
xmax = max(abs(xd)); % this takes into account penalty due to enhanced PAPR after pulse shapping 
xd = Tx.Mod.Vgain*Tx.Mod.Vswing/2*xd/xmax + Tx.Mod.Vbias; % Normalized to have excursion from approx +-Vswing 
% Note: Vswing and Vswing and Vbias are normalized by Vpi/2. This scaling
% and offset is done before the DAC in order to have right pre-distortion
% function, and not to have the right voltage values that will drive the
% modulator

%% Predistortion to compensate for MZM non-linear response
if isfield(sim, 'mzm_predistortion') && strcmpi(sim.mzm_predistortion, 'analog')    
    xd = 2/pi*asin(sqrt(abs(xd))).*sign(xd); % apply predistortion
end  

%% ================================ DAC ===================================
% Set DAC time offset in order to remove group delay due to pulse shaping. 
% This way the first sample of xt will be the center of the first pulse. 
% This is only important for plotting.
Tx.DAC.offset = sim.Mct/mpam.pulse_shape.sps*(length(mpam.pulse_shape.h)-1)/2;
xt = dac(xd, Tx.DAC, sim, sim.shouldPlot('DAC output'));
% Note: Driving signal xd must be normalized by Vpi

%% ============================= Modulator ================================
Tx.Laser.PdBm = Watt2dBm(Tx.Ptx);
Ecw = Tx.Laser.cw(sim);
Etx = mzm(Ecw, xt, Tx.Mod); % transmitted electric field

% Adjust power to make sure desired power is transmitted
Etx = Etx*sqrt(Tx.Ptx/mean(abs(Etx).^2));

P = abs(Etx).^2;
rexdB = 10*log10(min(P(1:sim.Mct:end))/max(P(1:sim.Mct:end)))

%% ========================= Fiber propagation ============================
Erx = Etx;
% link_gain = Amp.Gain*Rx.PD.R;
for k = 1:length(Fibers)
    fiberk = Fibers(k); 
    
    Erx = fiberk.linear_propagation(Erx, sim.f, Tx.Laser.wavelength); % propagation through kth fiber in Fibers
end

%% ========================= Preamplifier =================================
if isfield(sim, 'preAmp') && sim.preAmp
    [Erx, OSNRdBest] = Amp.amp(Erx, sim.fs);
    
    %% Optical bandpass filter
%     Hopt = ifftshift(Rx.optfilt.H(sim.f/sim.fs));
%     Erx = [ifft(fft(Erx(1, :)).*Hopt);...
%         ifft(fft(Erx(2, :)).*Hopt)];
else
    Erx = [Erx; zeros(size(Erx))];
end

% Ensure that signal has the desired output power (EDFA in power mode)
Att = dBm2Watt(sim.PrxdBm)/dBm2Watt(power_meter(Erx));
Erx = Erx*sqrt(Att);  % keep constant received power of 

% Measure OSNR
OSNRdB = estimate_osnr(Erx, Tx.Laser.wavelength, sim.f, sim.shouldPlot('OSNR'))


%% ========================== Receiver ====================================

if isfield(sim, 'coherentReceiver') && sim.coherentReceiver
    %% Do direct detection using coherent receiver
    yk = coherent_imdd_rx(Erx, mpam, Tx.Laser, Rx, sim); 
    filt = design_filter('butter', 5, 0.7*mpam.Rs/(sim.fs/2));
    Hrx = filt.H(sim.f/sim.fs);
else
    %% Direct detection
    % Direct detection and add thermal noise
    % PD.detect(Ein: Input electric field, fs: sampling rate of samples in Ein,
    % noise statistics {'gaussian', 'no noise'}, N0: one-sided PSD of thermal
    % noise)
    yt = Rx.PD.detect(Erx(1, :), sim.fs, 'gaussian', Rx.N0);

    % Gain control
    yt = yt - mean(yt);
    yt = yt*sqrt(mean(abs(mpam.a).^2)/(sqrt(2)*mean(abs(yt).^2)));
    
    %% ADC
    % ADC performs filtering, quantization, and downsampling
    % For an ideal ADC, ADC.ENOB = Inf
    % Rx.ADC.offset = -1;       
    Hrx = Rx.ADC.filt.H(sim.f/sim.fs);
    [yk, ~, ytf] = adc(yt, Rx.ADC, sim);
end

%% =========================== Equalization ===============================
Rx.eq.trainSeq = dataTX; % training sequence
[yd, Rx.eq] = equalize(Rx.eq, yk, [], mpam, sim, sim.shouldPlot('Equalizer'));

% Symbols to be discard in BER calculation
ndiscard = [1:Rx.eq.Ndiscard(1)+sim.Ndiscard (sim.Nsymb-Rx.eq.Ndiscard(2)-sim.Ndiscard):sim.Nsymb];
ydfull = yd;
yd(ndiscard) = []; 
dataTX(ndiscard) = [];

%% =========================== Demodulation ===============================
% dataRX = mpam.demod(yd);
[dataRX, mpam] = mpam.demod_sweeping_thresholds(yd, dataTX);

%% Counted BER
[~, ber_count] = biterr(dataRX, dataTX)

%% ====================== Gaussian approxiation ===========================
% Note: this calculation includes noise enhancement due to equalization,
% but dominant noise in pre-amplified system is the signal-spontaneous beat
% noise, which is not Gaussian distributed

% Noise bandwidth
Hrxeq = abs(Hrx.*Rx.eq.Hff(sim.f/(Rx.eq.ros*mpam.Rs))).^2;
Hrxeq = Hrxeq/interp1(sim.f, Hrxeq, 0); % normalize to have unit gain at DC
noiseBW = 0.5*trapz(sim.f, Hrxeq);

BWopt = Rx.optfilt.noisebw(sim.fs); % optical filter noise bandwidth; Not divided by 2 because optical filter is a bandpass filter

% Noise std for intensity level Plevel
Npol = 2; % number of polarizations. Npol = 1, if polarizer is present, Npol = 2 otherwise.
if isfield(sim, 'coherentReceiver') && sim.coherentReceiver % includes shot noise penalty
    noise_std = @(Plevel) sqrt(2*Att*Plevel*Amp.N0*noiseBW + 1/2*Npol*Amp.N0^2*BWopt*noiseBW...
        + 8*Plevel*Rx.PD.varShot(Rx.PlodBm, noiseBW));
else
    noise_std = @(Plevel) sqrt(2*Att*Plevel*Amp.N0*noiseBW + 1/2*Npol*Amp.N0^2*BWopt*noiseBW);
end
    
% Note: Plevel is divided by amplifier gain to obtain power at the amplifier input

% AWGN approximation
Prx = dBm2Watt(sim.PrxdBm) - Amp.N0*sim.fs;
% Note: fs is not divided by 2 to account for noise in two pols
mpamRef = mpamRef.adjust_levels(Prx, rexdB);

ber_awgn = mpamRef.berAWGN(noise_std);

%% 4-PAM Bias analysis
Vset = [];
if mpam.M == 4 && strcmpi(sim.mzm_predistortion, 'none')
	p(1) = mean(yd(yd < mpam.b(1)));
    p(2) = mean(yd(yd > mpam.b(1) & yd < mpam.b(2)));
    p(3) = mean(yd(yd > mpam.b(2) & yd < mpam.b(3)));
    p(4) = mean(yd(yd > mpam.b(3)));
    
    mzm_nonlinearity = @(levels, V) V(3)*abs(sin(pi/2*(levels*V(1) + V(2)))).^2 + V(4);
    
    [Vset, fval, exitflag] = fminsearch(@(V) norm(p.' - (mzm_nonlinearity(mpam.a, V) - mean(mzm_nonlinearity(mpam.a, V)))), [0.5 0.5 1 0]);
    
    if exitflag ~= 1
        disp('4-PAM bias control did not converge')
    end
    
    fprintf('Vgain/Vpi = %.2f\n', Vset(1))
    fprintf('Vbias/Vpi = %.2f\n', Vset(2))
end

%% Plots
if sim.shouldPlot('Optical eye diagram')  
    Ntraces = 500;
    Nstart = sim.Ndiscard*sim.Mct + 1;
    Nend = min(Nstart + Ntraces*2*sim.Mct, length(Etx));
    figure(103), clf
    subplot(121), box on
    eyediagram(abs(Etx(Nstart:Nend)).^2, 2*sim.Mct)
    title('Transmitted optical signal eye diagram')
    subplot(122), box on
    eyediagram(abs(Erx(1, Nstart:Nend)).^2, 2*sim.Mct)
    title('Received optical signal eye diagram')
    drawnow
end

if sim.shouldPlot('Received signal eye diagram')
    mpam = mpam.norm_levels();
    Ntraces = 500;
    Nstart = sim.Ndiscard*sim.Mct + 1;
    Nend = min(Nstart + Ntraces*2*sim.Mct, length(ytf));
    figure(104), clf, box on, hold on
    eyediagram(ytf(Nstart:Nend), 2*sim.Mct)
    title('Received signal eye diagram')
    a = axis;
    h1 = plot(a(1:2), mpam.a*[1 1], '-k');
    h2 = plot(a(1:2), mpam.b*[1 1], '--k');
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
