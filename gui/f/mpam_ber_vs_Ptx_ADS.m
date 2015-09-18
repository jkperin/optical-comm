function ber_ads = mpam_ber_vs_Ptx_ADS(mpam, tx, fiber1, soa1, apd1, rx, sim)

assert(mpam.M == 4, 'ADS co-simulation only supports PAM 4.')

dBm2Watt = @(x) 1e-3*10.^(x/10);

% Load file with ADS simulation results
try
    load(sim.ads.filename)
    
    % Simulation parameters
    sim.Nsymb = Nsymb; % Number of symbols in montecarlo simulation
    sim.Mct = Mct;     % Oversampling ratio to simulate continuous time (must be odd) 
    sim.N = Nsymb*sim.Mct; % number points in 'continuous-time' simulation
    
    % If they don't exist an error will be thrown
    Pout;
    Poutnf;    
catch e
    disp('Invalid ADS file!')
    rethrow(e);
end

% Optimize levels and decision instant
[Pthresh, Plevels, OptDecisionInstant] = PAM4_timing_opt(Pout, Mct, false);

% Set optimized levels
mpam = mpam.set_levels(Plevels, Pthresh);

% Shift signal so that optimal decision threshold appears in the middle of
% the pulse (required for matched filtering)
Poutnf = circshift(Poutnf, [0, (Mct+1)/2 - OptDecisionInstant-1]);

%% Matlab Simulation
% Time and frequency
sim.fs = mpam.Rs*Mct;  % sampling frequency in 'continuous-time'

dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

sim.t = t;
sim.f = f;

% Transmitter
tx.rexdB = 10*log10(mpam.a(1)/mpam.a(end));          % exctinction ratio in dB. Defined as Pmin/Pmax

tx = rmfield(tx, 'modulator');
tx_ads = tx;

% Modulator frequency response
tx.modulator.fc = 24e9; % modulator cut off frequency
tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
tx.modulator.h = @(t) [0*t(t < 0) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0))];
tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

% Receiver
if isempty(apd1) % use pin photodetector
    photodiode = apd(0, 0, Inf, rx.R, rx.Id);
else
    photodiode = apd1;
end
    
if isempty(soa1)
    Gsoa = 1;
else
    Gsoa = soa1.Gain;
end

%% Receiver
if strcmp(rx.filter.type, 'matched')
    % Electric Lowpass Filter
    rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);
else
    rx.elefilt = design_filter(rx.filter.type, rx.filterN, rx.filterBw/(sim.fs/2));
end

filt = design_filter('bessel', 5, 19e9/(sim.fs/2));

%% Channel Response
% Hch does not include transmitter or receiver filter
Hch = tx.modulator.H(sim.f).*exp(1j*2*pi*sim.f*tx.modulator.grpdelay)...
    .*fiber1.H(sim.f, tx).*photodiode.H(sim.f);

link_gain = Gsoa*photodiode.R*photodiode.Gain*fiber1.link_attenuation(tx.lamb); % Overall link gain. Equivalent to Hch(0)

% Symbols to be discard in BER calculation
Ndiscard = 50*Mct*[1 1]; % discard symbols to account for transients in ADS simulation
if rx.eq.adaptive
    Ndiscard(1) = Ndiscard(1) + rx.eq.Ntrain;
end
if isfield(rx.eq, 'Ntaps')
    Ndiscard = Ndiscard + rx.eq.Ntaps;
end
ndiscard = [1:Ndiscard(1) sim.Nsymb-Ndiscard(2):sim.Nsymb];

% ndiscard = [1:(sim.Ndiscard+eq.Ntrain) (sim.Nsymb-sim.Ndiscard-eq.Ntrain):sim.Nsymb];

% Equalization
rx.eq.TrainSeq = dataTXref;

Ptx = dBm2Watt(tx.PtxdBm);

ber_ads = zeros(sim.ads.Nrealizations, length(Ptx));
for kk = 1:sim.ads.Nrealizations
    for k = 1:length(Ptx)
        tx.Ptx = Ptx(k);

        % Ajust levels to desired transmitted power and extinction ratio
        Ptads = Poutnf.'/mean(Poutnf)*tx.Ptx;
        Pmax = Plevels(end)/mean(Poutnf)*tx.Ptx;
        mpam = mpam.adjust_levels(tx.Ptx, tx.rexdB);
        Pthresh = mpam.b;
        
        Ptads([sim.Mct*Ndiscard(1) end-sim.Mct*Ndiscard(2):end]) = 0; % zero sim.Ndiscard last symbbols

        % Generate optical signal
        [Etads, ~] = optical_modulator(Ptads, tx_ads, sim); % without frequency response, which is already included in ADS simulation

        % Fiber propagation
        [Etads, Ptads] = fiber1.linear_propagation(Etads, sim.f, tx.lamb);
        
        if ~isempty(soa1)
            % Amplifier
            et = soa1.amp(Etads, sim.fs);

            % Optical bandpass filter
            eo = [ifft(fft(et(:, 1)).*ifftshift(rx.optfilt.H(f))),...
                ifft(fft(et(:, 2)).*ifftshift(rx.optfilt.H(f)))];

            % Direct detection and add noises
            if isfield(sim, 'polarizer') && ~sim.polarizer 
                Ptads = abs(eo(:, 1)).^2 + abs(eo(:, 2)).^2;
            else % by default assumes that polarizer is used
                Ptads = abs(eo(:, 1)).^2;
            end
        end

        % Detect and add noises
        ytads = photodiode.detect(Ptads, sim.fs, 'gaussian', rx.N0);

        if ~strcmp(rx.eq.type, 'None')
            ytads = real(ifft(fft(ytads).*ifftshift(filt.H(sim.f/sim.fs))));
        end

        % Automatic gain control
%         Pmax = 2*tx.Ptx/(1 + 10^(-abs(tx.rexdB)/10));
        ytads = ytads/(Pmax*link_gain);
        mpam = mpam.norm_levels;

        % Equalize       
        ydads = equalize(rx.eq, ytads, Hch, mpam, rx, sim);

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
        [~, ber_ads(kk, k)] = biterr(dataRXads, dataTX);
    end
end

ber_ads = mean(ber_ads, 1);

if sim.ads.eyediagram
    n = Ndiscard(1):length(t)-Ndiscard(2);
    n = n(1:min(2^14, length(n)));
        
    eyediagram(circshift(ytads(n), [-(Mct+1)/2+1, 0]), 2*sim.Mct)
    xlabel('T/(2T_s)')
    ylabel('Received Electric Signal')
    axis([-0.5, 0.5 0 1.2])
    hold on
    plot([-0.5 0.5], mpam.b*[1 1], 'k', 'LineWidth', 2)
    title('Received Electric Signal')
end