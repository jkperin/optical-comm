function ber_ads = mpam_ber_vs_Ptx_ADS(tx, fiber1, rx, sim)
dBm2Watt = @(x) 1e-3*10.^(x/10);

load SampleData4PAMlong

%% Matlab Simulation
% Simulation parameters
sim.Nsymb = Nsymb; % Number of symbols in montecarlo simulation
sim.Mct = Mct_odd;     % Oversampling ratio to simulate continuous time (must be odd) 
sim.N = Nsymb*sim.Mct; % number points in 'continuous-time' simulation

% Time and frequency
sim.fs = mpam.Rs*Mct;  % sampling frequency in 'continuous-time'

dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

sim.t = t;
sim.f = f;

% Transmitter
tx.rexdB = -5;          % exctinction ratio in dB. Defined as Pmin/Pmax

% Modulator frequency response
tx.modulator.fc = 24e9; % modulator cut off frequency
tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
tx.modulator.h = @(t) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0));
tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

tx_ads = tx;
tx_ads = rmfield(tx_ads, 'modulator');

% Receiver
pin = apd(0, 0, Inf, rx.R, rx.Id);

%% Receiver
if strcmp(rx.filter.type, 'matched')
    % Electric Lowpass Filter
    rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);
else
    rx.elefilt = design_filter(rx.filter.type, rx.filterN, rx.filterBw/(sim.fs/2));
end

filt = design_filter('bessel', 5, 19e9/(sim.fs/2));

% Equalization
eq = rx.eq;

% Overall link gain
link_gain = fiber1.link_attenuation(tx.lamb)*pin.R;

% Symbols to be discard in BER calculation
Ndiscard = sim.Ndiscard*[1 1];
if isfield(rx.eq, 'Ntrain')
    Ndiscard(1) = Ndiscard(1) + rx.eq.Ntrain;
end
if isfield(rx.eq, 'Ntaps')
    Ndiscard = Ndiscard + rx.eq.Ntaps;
end
ndiscard = [1:Ndiscard(1) sim.Nsymb-Ndiscard(2):sim.Nsymb];

% ndiscard = [1:(sim.Ndiscard+eq.Ntrain) (sim.Nsymb-sim.Ndiscard-eq.Ntrain):sim.Nsymb];

%
% Equalization
rx.eq.TrainSeq = dataTXref;

Poutnf_odd = [0 Poutnf_odd];

Ptx = dBm2Watt(tx.PtxdBm);
for k = 1:length(Ptx)
    tx.Ptx = Ptx(k);

	% Ajust levels to desired transmitted power and extinction ratio
    Ptads = Poutnf_odd.'/mean(Poutnf_odd)*tx.Ptx;
    mpam.adjust_levels(tx.Ptx, tx.rexdB);
    
    Ptads([sim.Mct*Ndiscard(1) end-sim.Mct*Ndiscard(2):end]) = 0; % zero sim.Ndiscard last symbbols

	% Generate optical signal
    [Etads, ~] = optical_modulator(Ptads, tx_ads, sim); % without frequency response, which is already included in ADS simulation

	% Fiber propagation
    [~, Ptads] = fiber1.linear_propagation(Etads, sim.f, tx.lamb);

    % Detect and add noises
    ytads = pin.detect(Ptads, sim.fs, 'gaussian', rx.N0);
    
    if ~strcmp(rx.eq.type, 'None')
        ytads = real(ifft(fft(ytads).*ifftshift(filt.H(sim.f/sim.fs))));
    end
  
    % Automatic gain control
    Pmax = 2*tx.Ptx/(1 + 10^(-abs(tx.rexdB)/10));
    ytads = ytads/(Pmax*link_gain);
    mpam.norm_levels;
    
    % Equalization
    ydads = equalize(rx.eq.type, ytads, mpam, tx, fiber1, rx, sim);

    ydads(ndiscard) = [];
    dataTX = dataTXref;
    dataTX(ndiscard) = [];
    
    % Detect ADS signal
    if mpam.M == 2
        ydads = ydads - mean(ytads);
        dataRXads = zeros(size(ydads)).';
        dataRXads(ydads > 0) = 1;
    else        
%         mpam.b = Pthresh.'/mean(Poutnf)*tx.Ptx;
%         dataRXads = mpam.demod(ydads);
%         mpam.b = btemp;
        dataRXads = mpam.demod(ydads);
    end
    
    % BER
    [~, ber_ads(k)] = biterr(dataRXads, dataTX);
end

if sim.ads.eyediagram
    n = Ndiscard(1):length(t)-Ndiscard(2);
    n = n(1:min(2^14, length(n)));
    
    eyediagram(abs(Etads(n)).^2, 2*sim.Mct)
    xlabel('T/(2T_s)')
    ylabel('Transmitted Optical Signal')
    axis([-0.5, 0.5 0 1.2*Pmax])
    title(sprintf('Transmitted Optical Signal: Ptx = %.2f mW, RIN = %.1f dB/Hz, ER = %.1f dB', tx.Ptx*1e3, tx.RIN, tx.rexdB))
    
    eyediagram(ytads(n), 2*sim.Mct)
    xlabel('T/(2T_s)')
    ylabel('Received Electric Signal')
    axis([-0.5, 0.5 0 1.2])
    hold on
    plot([-0.5 0.5], mpam.b*[1 1], 'k')
    title('Received Electric Signal')
end