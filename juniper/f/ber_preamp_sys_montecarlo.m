function [ber_count, ber_awgn] = ber_preamp_sys_montecarlo(mpam, tx, Fibers, Amp, rx, sim)
%% Calculate BER of pre-amplified IM-DD system through montecarlo simulation
% Inputs:
% - mpam: PAM class
% - tx: struct with transmitter paramters
% - Fibers: array of Fiber classes containing the fibers used in
% transmittion i.e., Fibers = [SMF DCF] or Fibers = [SMF];
% - Amp: pre-amplifier using SOA class
% - rx: struct with receiver parameters
% - sim: struct with simulation parameters

% Normalized frequency
f = sim.f/sim.fs;

% Modulator frequency response
if isfield(tx, 'modulator')
    Hmod = tx.modulator.H(sim.f).*exp(1j*2*pi*sim.f*tx.modulator.grpdelay);
else
    Hmod = 1;
end

% Photodiode frequency response
Hpd = rx.PD.H(sim.f);

% Ajust levels to desired transmitted power and extinction ratio
mpam = mpam.adjust_levels(tx.Ptx, tx.rexdB);
Pmax = mpam.a(end); % used in the automatic gain control stage

% Modulated PAM signal
dataTX = randi([0 mpam.M-1], 1, sim.Nsymb); % Random sequence
xt = mpam.mod(dataTX, sim.Mct).';
xt(1:sim.Mct*sim.Ndiscard) = 0; % zero sim.Ndiscard first symbols
xt(end-sim.Mct*sim.Ndiscard+1:end) = 0; % zero sim.Ndiscard last symbbols

% Generate optical signal
[Et, ~] = optical_modulator(xt, tx, sim);

% Fiber propagation
Hch = Hmod.*Hpd;
link_gain = Amp.Gain*rx.PD.R;
for k = 1:length(Fibers)
    fiberk = Fibers(k);
    
    link_gain = link_gain*fiberk.link_attenuation(tx.lamb); % Overall link gain
    
    Hch = Hch.*fiberk.H(sim.f, tx); % frequency response of the channel (used in designing the equalizer)
    
    Et = fiberk.linear_propagation(Et, sim.f, tx.lamb); % propagation through kth fiber in Fibers
end

% Amplifier
Et = Amp.amp(Et, sim.fs);

% Optical bandpass filter
Eo = [ifft(fft(Et(1, :)).*ifftshift(rx.optfilt.H(f)));...
    ifft(fft(Et(2, :)).*ifftshift(rx.optfilt.H(f)))];

%% Direct detection and add thermal noise
% PD.detect(Ein: Input electric field, fs: sampling rate of samples in Ein,
% noise statistics {'gaussian', 'no noise'}, N0: one-sided PSD of thermal
% noise)
yt = rx.PD.detect(Eo, sim.fs, 'gaussian', rx.N0);

% Automatic gain control
% Pmax = 2*tx.Ptx/(1 + 10^(-abs(tx.rexdB)/10)); % calculated from mpam.a
yt = yt/(Pmax*link_gain); % just refer power values back to transmitter
mpam = mpam.norm_levels;

%% Equalization
if isfield(rx, 'eq') && isfield(tx, 'modulator')
    rx.eq.TrainSeq = dataTX;
else % otherwise only filter using rx.elefilt
    rx.eq.type = 'None';
end
   
% Equalize
[yd, rx.eq] = equalize(rx.eq, yt, Hch, mpam, rx, sim);

% Symbols to be discard in BER calculation
ndiscard = [1:rx.eq.Ndiscard(1) sim.Nsymb-rx.eq.Ndiscard(2):sim.Nsymb];
yd(ndiscard) = []; 
dataTX(ndiscard) = [];

% Demodulate
dataRX = mpam.demod(yd.');

% Counted BER
[~, ber_count] = biterr(dataRX, dataTX);

%% AWGN approximation
% Note: this calculation includes noise enhancement due to equalization,
% but dominant noise in pre-amplified system is the signal-spontaneous beat
% noise, which is not Gaussian distributed

% Noise bandwidth
if isfield(rx, 'eq')
    noiseBW = trapz(sim.f, abs(rx.eq.Hrx.'.*rx.eq.Hff(sim.f/(mpam.Rs*rx.eq.ros))).^2)/2;
else 
    noiseBW = rx.elefilt.noisebw(sim.fs)/2;
end
fprintf('Noise bandwidth = %.2f GHz\n', noiseBW/1e9);

BWopt = rx.optfilt.noisebw(sim.fs); % optical filter noise bandwidth; Not divided by 2 because optical filter is a bandpass filter

% Thermal noise
varTherm = rx.N0*noiseBW; % variance of thermal noise

% RIN
if isfield(sim, 'RIN') && sim.RIN
    varRIN =  @(Plevel) 10^(tx.RIN/10)*Plevel.^2*noiseBW;
else
    varRIN = @(Plevel) 0;
end

% Noise std for intensity level Plevel
Npol = 2; % number of polarizations. Npol = 1, if polarizer is present, Npol = 2 otherwise.
noise_std = @(Plevel) sqrt(varTherm + rx.PD.varShot(Plevel, noiseBW) + rx.PD.R^2*varRIN(Plevel)...
    + rx.PD.R^2*Amp.var_awgn(Plevel/Amp.Gain, noiseBW, BWopt, Npol));
% Note: Plevel is divided by amplifier gain to obtain power at the amplifier input

% AWGN approximation
mpam = mpam.adjust_levels(tx.Ptx*link_gain, tx.rexdB);

ber_awgn = mpam.ber_awgn(noise_std);

%% Plots
if isfield(sim, 'Plots') && sim.Plots.isKey('Equalizer')  && sim.Plots('Equalizer')
    figure(100)
    for k = 1:size(rx.eq.num, 2)
        [h, w] = freqz(rx.eq.num(:, k), rx.eq.den(:, k));
        subplot(121), hold on, box on
        plot(w/(2*pi), abs(h).^2)
        
        subplot(122), hold on, box on
        plot(w/(2*pi), unwrap(angle(h)))
        
        subplot(121), hold on, box on
        xlabel('Normalized frequency')
        ylabel('|H(f)|^2')
        
        subplot(122), hold on, box on
        xlabel('Normalized frequency')
        ylabel('arg(H(f))')
    end
    drawnow
end

if isfield(sim, 'Plots') && sim.Plots.isKey('Output signal')  && sim.Plots('Output signal')   
    figure(101), hold on
    plot(link_gain*Pt)
    plot(yt, '-k')
    plot(ix, yd, 'o')
    legend('Transmitted power', 'Received signal', 'Samples')
    drawnow
end

if isfield(sim, 'Plots') && sim.Plots.isKey('Heuristic noise pdf') && sim.Plots('Heuristic noise pdf')
    figure(102)
    [nn, xx] = hist(yd(dataTX == 2), 50);
    nn = nn/trapz(xx, nn);
    bar(xx, nn)
    drawnow
end
