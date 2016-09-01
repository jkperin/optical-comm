function [ber_count, ber_awgn, OSNRdB] = ber_preamp_sys_montecarlo(mpam, Tx, Fibers, Amp, Rx, sim)
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
    
    xd = mpamdb.mod(dataTX); % Modulated PAM signal

    xd_enc = duobinary_encoding(xd);
      
    predist = @(p) 2/pi*asin(sqrt(p));
    dist = @(v) sin(pi/2*v).^2;
    
    % 4 intensity levels are trnasmitted. Modulator bias must be set to 0
    Pswing = dist(Tx.Mod.Vswing/2);
    DP = Pswing/(mpam.M-1);
    Pk = 0:DP:Pswing;
    Vkp = predist(Pk);
    Vk = [-Vkp(end:-1:2) Vkp];
    Vk = Vk/Vk(end);
        
    ximp = upsample(Vk(xd_enc+mpam.M), mpam.pulse_shape.sps);
    xd = filter(ones(1, mpam.pulse_shape.sps)/mpam.pulse_shape.sps, 1, ximp);
elseif isfield(sim, 'mzm_predistortion') && strcmpi(sim.mzm_predistortion, 'levels') %% Ordinary PAM with predistorted levels
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
if all(isreal(xd))
    xt = dac(xd, Tx.DAC, sim, sim.shouldPlot('DAC output'));
else
    xti = dac(real(xd), Tx.DAC, sim, sim.shouldPlot('DAC output'));
    xtq = dac(imag(xd), Tx.DAC, sim, sim.shouldPlot('DAC output'));
    xt = xti + 1j*xtq;
end

%% Driver
% Adjust gain to compensate for preemphasis
xt = Tx.Vgain*(xt - mean(xt)) + Tx.VbiasAdj*mean(xt);

%% ============================= Modulator ================================
Tx.Laser.PdBm = Watt2dBm(Tx.Ptx);
Ecw = Tx.Laser.cw(sim);
Etx = mzm(Ecw, xt, Tx.Mod); % transmitted electric field

% Adjust power to make sure desired power is transmitted
Etx = Etx*sqrt(Tx.Ptx/mean(abs(Etx).^2));

% Chirp
if isfield(Tx.Mod, 'alpha') && Tx.Mod.alpha ~= 0
    disp('Chirp added')
    Etx = Etx.*exp(1j*Tx.Mod.alpha/2*log(abs(Etx).^2));
end

% Calculate extinction ratio and intensity levels position
P = abs(Etx).^2;
Psamp = P(1:sim.Mct:end);
for k = 1:mpam.M
    Pl(k) = mean(Psamp(dataTX == k-1));
end
Pl = sort(Pl, 'ascend');
rexdB = 10*log10(min(Pl)/max(Pl));
fprintf('> Estimated extinction ratio = %.2f dB\n', rexdB)

%% ========================= Fiber propagation ============================
Erx = Etx;
% link_gain = Amp.Gain*Rx.PD.R;
for k = 1:length(Fibers)
    fiberk = Fibers(k); 
    
    Erx = fiberk.linear_propagation(Erx, sim.f, Tx.Laser.wavelength); % propagation through kth fiber in Fibers
end

%% ========================= Preamplifier =================================
if isfield(sim, 'preAmp') && sim.preAmp
    [Erx, OSNRdB.theory] = Amp.amp(Erx, sim.fs);
    
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
OSNRdB.measured = estimate_osnr(Erx, Tx.Laser.wavelength, sim.f, sim.shouldPlot('OSNR'))

%% ========================== Receiver ====================================
if isfield(sim, 'coherentReceiver') && sim.coherentReceiver
    %% Do direct detection using coherent receiver
    disp('!! Using coherent receiver to do direct detection');
    filt = design_filter('butter', 5, mpam.Rs/(sim.fs/2));
    Hrx = filt.H(sim.f/sim.fs);
    Rx.Filt = filt;
    yk = coherent_imdd_rx(Erx, mpam, Tx.Laser, Rx, sim); 
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
%     [yd, Rx.eq] = nonlinear_equalizer(Rx.eq, yk, mpam, sim, sim.shouldPlot('Equalizer'));

% Symbols to be discard in BER calculation
ndiscard = [1:Rx.eq.Ndiscard(1)+sim.Ndiscard (sim.Nsymb-Rx.eq.Ndiscard(2)-sim.Ndiscard):sim.Nsymb];
ydfull = yd;
yd(ndiscard) = []; 
dataTX(ndiscard) = [];

%% =========================== Detection ==============================
% dataRX = mpam.demod(yd);
[dataRX, mpam] = mpam.demod_sweeping_thresholds(yd, dataTX);

%% Counted BER
[~, ber_count] = biterr(dataRX, dataTX);

%% ====================== Gaussian approximation ===========================
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
    noise_std = @(Plevel) sqrt(Amp.varSigSpont(Att*Plevel/Amp.Gain, noiseBW)... % sig-spont
        + Amp.varSpontSpont(noiseBW, BWopt, Npol)... % spont-spont
        + 8*Plevel*Rx.PD.varShot(dBm2Watt(Rx.PlodBm), noiseBW)); % LO shot noise
else
    noise_std = @(Plevel) sqrt(Amp.varSigSpont(Att*Plevel/Amp.Gain, noiseBW)... % sig-spont
        + Amp.varSpontSpont(noiseBW, BWopt, Npol)); % spont-spont
end
    
% Note: Plevel is divided by amplifier gain to obtain power at the amplifier input

% AWGN approximation
Prx = dBm2Watt(sim.PrxdBm) - Npol*Amp.Ssp*2*sim.fs/2; % Ssp is one-sided PSD per real dimension
mpamRef = mpamRef.adjust_levels(Prx, rexdB);

ber_awgn = mpamRef.berAWGN(noise_std);

%% 4-PAM Bias analysis
Vset = [];
if mpam.M == 4
    if isfield(sim, 'duobinary') && sim.duobinary
        disp('Duobinary')
        a = Vk(4:end).'/max(Vk);
    elseif ~strcmpi(sim.mzm_predistortion, 'none')
        disp('!! Levels were predistorted')
        mpamPredist = mpamPredist.unbias();
        a = mpamPredist.a;
    else
        a = mpam.a;
    end
            
	p(1) = mean(yd(yd < mpam.b(1)));
    p(2) = mean(yd(yd > mpam.b(1) & yd < mpam.b(2)));
    p(3) = mean(yd(yd > mpam.b(2) & yd < mpam.b(3)));
    p(4) = mean(yd(yd > mpam.b(3)));
    
    mzm_nonlinearity = @(levels, V) V(3)*abs(sin(pi/2*(levels*V(1) + V(2)))).^2 + V(4);
    
    [Vset, fval, exitflag] = fminsearch(@(V) norm(p.' - (mzm_nonlinearity(a, V) - mean(mzm_nonlinearity(a, V)))), [0.5 0.5 1 0]);
    
    if exitflag ~= 1
        disp('4-PAM bias control did not converge')
    end
    
    fprintf('Vswing/Vpi = %.2f | expected %.2f\n', 2*Vset(1), Tx.Mod.Vswing)
    fprintf('Vbias/Vpi = %.2f | expected %.2f\n', Vset(2), Tx.Mod.Vbias)

    Vdrive = a*Vset(1) + Vset(2);
    Pset = abs(sin(pi/2*Vdrive)).^2;
    
    if sim.shouldPlot('Bias analysis')
        figure(233), clf, hold on, box on
        t = linspace(0, 1);
        plot(t, sin(pi/2*t).^2, 'k');
        plot((Vdrive*[1 1]).', [zeros(1, mpam.M); Pset.'], 'k');
        plot([zeros(1, mpam.M); Vdrive.'], ([1; 1]*Pset.'), 'k');
        xlabel('Driving signal')
        ylabel('Resulting power levels')
        axis([0 1 0 1])
        drawnow
    end
end

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
