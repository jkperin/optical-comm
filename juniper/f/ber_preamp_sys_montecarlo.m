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

mpamRef = mpam; % used in Gaussian approximation

% Normalized frequency
f = sim.f/sim.fs;

dataTX = randi([0 mpam.M-1], 1, sim.Nsymb); % Random sequence   
Nzero = 10;
dataTX([1:Nzero end-Nzero+1:end]) = 0; % set first and last Nzero symbols to 0 to make sequence periodic

%% Duobinary enconding
if isfield(sim, 'duobinary') && sim.duobinary
    mpamdb = mpam.set_levels(0:mpam.M-1, 0.5 + (0:mpam.M-2));    
    
    xd = mpamdb.signal(dataTX); % Modulated PAM signal
    
    xd = duobinary_encoding(xd);
    
    xd_enc = xd;
else
    xd = mpam.signal(dataTX); % Modulated PAM signal
end  

%% Predistortion to compensate for MZM non-linear response
if isfield(sim, 'predistortion') && strcmpi(sim.predistortion, 'levels')    
    xd = sqrt(abs(xd)).*sign(xd); % apply predistortion
end  

%% Preemphasis
if sim.preemphasis
    femph = abs(freq_time(sim.Nsymb*sim.ros.txDSP, mpam.Rs*sim.ros.txDSP));
    femph(femph >= sim.preemphRange) = 0;
    emphasis_filter = 10.^(polyval([-0.0013 0.5846 1.5859], femph/1e9)/20);    

    xd = real(ifft(fft(xd).*ifftshift(emphasis_filter)));
end

%% Driver
xd = xd - mean(xd); % Remove DC bias, since modulator is IM-MZM
xmax = max(abs(xd)); % this takes into account penalty due to enhanced PAPR after pulse shapping 
xd = Tx.Mod.Vswing/2*xd/xmax + Tx.Mod.Vbias; % Normalized to have excursion from approx +-Vswing 
% Note: Vswing and Vswing and Vbias are normalized by Vpi/2. This scaling
% and offset is done before the DAC in order to have right pre-distortion
% function, and not to have the right voltage values that will drive the
% modulator

%% Predistortion to compensate for MZM non-linear response
if isfield(sim, 'predistortion') && strcmpi(sim.predistortion, 'analog')    
    xd = 2/pi*asin(sqrt(abs(xd))).*sign(xd); % apply predistortion
end  

%% DAC
% Set DAC time offset in order to remove group delay due to pulse shaping. 
% This way the first sample of xt will be the center of the first pulse. 
% This is only important for plotting.
Tx.DAC.offset = sim.Mct/mpam.pulse_shape.sps*(length(mpam.pulse_shape.h)-1)/2;

xt = dac(xd, Tx.DAC, sim);
% Note: Driving signal xd must be normalized by Vpi

%% Generate optical signal
Tx.Laser.PdBm = Watt2dBm(Tx.Ptx);
Ecw = Tx.Laser.cw(sim);
Etx = mzm(Ecw, xt, Tx.Mod);

% Chirp
% Adds transient chirp just to measure its effect. MZM in push pull should
% have no chirp
if isfield(Tx, 'alpha') && Tx.alpha ~= 0
    disp('chirp added!')
    Etx = Etx.*exp(1j*Tx.alpha/2*log(abs(Etx).^2));
end

% Adjust power to make sure desired power is transmitted
Etx = Etx*sqrt(Tx.Ptx/mean(abs(Etx).^2));

%% Fiber propagation
Erx = Etx;
Prx = Tx.Ptx;
for k = 1:length(Fibers)
    fiberk = Fibers(k); 
    
    Erx = fiberk.linear_propagation(Erx, sim.f, Tx.Laser.wavelength); % propagation through kth fiber in Fibers
    
    Prx = Prx*fiberk.link_attenuation(Tx.Laser.wavelength);
end

%% Pre-amplifier
[Erx, OSNRdB] = Amp.amp(Erx, sim.fs);
Prx = Prx*Amp.Gain;

%% Optical bandpass filter
Hopt = ifftshift(Rx.optfilt.H(f));
Erx = [ifft(fft(Erx(1, :)).*Hopt);...
    ifft(fft(Erx(2, :)).*Hopt)];

% Direct detection and add thermal noise
% PD.detect(Ein: Input electric field, fs: sampling rate of samples in Ein,
% noise statistics {'gaussian', 'no noise'}, N0: one-sided PSD of thermal
% noise)
yt = Rx.PD.detect(Erx, sim.fs, 'gaussian', Rx.N0);

%% Automatic gain control
mpam.b = mpam.b - mean(mpam.a);
mpam.a = mpam.a - mean(mpam.a);
mpam = mpam.norm_levels();

yt = yt - mean(yt);
yt = yt*sqrt(var(mpam.a)/(sqrt(2)*mean(abs(yt).^2)));

%% ADC
% ADC performs filtering, quantization, and downsampling
% For an ideal ADC, ADC.ENOB = Inf
% Rx.ADC.offset = -1;
switch lower(Rx.filtering)
    case 'antialiasing' % receiver filter is specified in ADC.filt
        Hrx = Rx.ADC.filt.H(f);
        [yk, ~, ytf] = adc(yt, Rx.ADC, sim);
    case 'matched' % receiver filter is matched filter
        Hrx = conj(Hch);
        [yk, ~, ytf] = adc(yt, Rx.ADC, sim, Hrx);
    otherwise
        error('ber_preamp_sys_montecarlo: Rx.filtering must be either antialiasing or matched')
end     

%% Equalization
Rx.eq.trainSeq = dataTX;
[yd, Rx.eq] = equalize(Rx.eq, yk, [], mpam, sim, sim.shouldPlot('Equalizer'));

% Symbols to be discard in BER calculation
ndiscard = [1:Rx.eq.Ndiscard(1)+sim.Ndiscard (sim.Nsymb-Rx.eq.Ndiscard(2)-sim.Ndiscard):sim.Nsymb];
ydfull = yd;
yd(ndiscard) = []; 
dataTX(ndiscard) = [];

%% Demodulate
% dataRX = mpam.demod(yd);
[dataRX, mpam] = mpam.demod_sweeping_thresholds(yd, dataTX);

%% Counted BER
[~, ber_count] = biterr(dataRX, dataTX);

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
noise_std = @(Plevel) sqrt(2*Plevel*Amp.N0*noiseBW + 1/2*Npol*Amp.N0^2*BWopt*noiseBW);
% Note: Plevel is divided by amplifier gain to obtain power at the amplifier input

% AWGN approximation
mpamRef = mpamRef.adjust_levels(Prx, -Inf);

ber_awgn = mpamRef.berAWGN(noise_std);

%% Plots
if sim.shouldPlot('Optical eye diagram')  
    Ntraces = 500;
    Nstart = sim.Ndiscard*sim.Mct + 1;
    Nend = min(Nstart + Ntraces*2*sim.Mct, length(Etx));
    figure(103), clf, box on
    eyediagram(abs(Etx(Nstart:Nend)).^2, 2*sim.Mct)
    title('Optical eye diagram')
    drawnow
end

if sim.shouldPlot('Received signal eye diagram')
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
    figure(105), clf, box on, hold on
    h1 = plot(ydfull, 'o');
    a = axis;
    h2= plot(a(1:2), (mpam.a*[1 1]).', '-k');
    h3 = plot(a(1:2), (mpam.b*[1 1]).', '--k');
    h4 = plot(Rx.eq.Ndiscard(1)*[1 1], a(3:4), ':k');
    h5 = plot((sim.Nsymb-Rx.eq.Ndiscard(2))*[1 1], a(3:4), ':k');
    legend([h1 h2(1) h3(1) h4], {'Equalized samples', 'PAM levels',...
        'Decision thresholds', 'BER measurement window'})
    title('Signal after equalization')
    axis([1 sim.Nsymb -1.2 1.2])
    drawnow
end

if sim.shouldPlot('Heuristic noise pdf')
    figure(106)
    [nn, xx] = hist(yd(dataTX == 2), 50);
    nn = nn/trapz(xx, nn);
    bar(xx, nn)
    drawnow
end
