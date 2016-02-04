%% Calculate BER for unamplified IM-DD link with APD
% bergauss = BER using gaussian approximation
% bertail_levels = BER of each level individually. This actually
% corresponds bertail_levels = p(error | given symbol)p(symbol)

function [bergauss, bergauss_levels, ber_awgn] = ber_apd_gauss(mpam, tx, fiber, apd, rx, sim)
verbose = sim.verbose;
sim.verbose = max(sim.verbose-1, 0);

Nsymb = mpam.M^sim.L; % number of data symbols 
% Number of symbols to be discared from begin and end of sequence
if isfield(rx, 'eq') && isfield(rx.eq, 'Ntaps')
    Ndisc = max(rx.eq.Ntaps, sim.L);
else
    Ndisc = sim.L;
end
N = sim.Mct*(Nsymb + 2*Ndisc); % total number of points 

% Frequency
df = 1/N;
f = (-0.5:df:0.5-df).';
sim.f = f*sim.fs; % redefine frequency to be used in optical_modulator.m

%% Channel Response
% Hch does not include transmitter or receiver filter
if isfield(tx, 'modulator')
    Hch = tx.modulator.H(sim.f).*exp(1j*2*pi*sim.f*tx.modulator.grpdelay)...
    .*fiber.H(sim.f, tx).*apd.H(sim.f);
else
    Hch = fiber.H(sim.f, tx).*apd.H(sim.f);
end

link_gain = apd.Gain*apd.R*fiber.link_attenuation(tx.lamb); % Overall link gain. Equivalent to Hch(0)

% Ajust levels to desired transmitted power and extinction ratio
mpam = mpam.adjust_levels(tx.Ptx, tx.rexdB);
Pmax = mpam.a(end); % used in the automatic gain control stage

% Modulated PAM signal
dataTX = debruijn_sequence(mpam.M, sim.L).'; % de Bruijin sequence
dataTXext = wextend('1D', 'ppd', dataTX, Ndisc); % periodic extension
xt = mpam.mod(dataTXext, sim.Mct);

% Generate optical signal
if ~isfield(sim, 'RIN')
    sim.RIN = false;
end

RIN = sim.RIN;
sim.RIN = false; % RIN is not modeled here since number of samples is not high enough to get accurate statistics
[Et, ~] = optical_modulator(xt, tx, sim);

% Fiber propagation
[~, Pt] = fiber.linear_propagation(Et, sim.f, tx.lamb);

% Direct detect
yt = apd.detect(Pt, sim.fs, 'no noise');

% Noise whitening filter
if sim.WhiteningFilter
    [Hw, yt] = apd.Hwhitening(sim.f, tx.Ptx, rx.N0, yt);
else
    Hw = 1;
end

%% Automatic gain control
% Normalize signal so that highest level is equal to 1
yt = yt/(Pmax*link_gain); 
mpam = mpam.norm_levels();

%% Equalization
[yd, rx.eq] = equalize(rx.eq, yt, Hw.*Hch, mpam, rx, sim);

% Symbols to be discard in BER calculation
yd = yd(Ndisc+1:end-Ndisc);

%% Detection
Pthresh = mpam.b; % decision thresholds referred to the receiver

% Noise std: includes RIN, shot and thermal noise (assumes gaussian stats)
noise_std = apd.stdNoise(Hw.*rx.eq.Hrx, rx.eq.Hff(sim.f/mpam.Rs), rx.N0, tx.RIN, sim);

%% Calculate signal-dependent noise variance after matched filtering and equalizer 
Ssh = apd.varShot(Pt, 1)/2; % two-sided shot noise PSD

% Receiver filter
% For symbol-rate sampling linear equalizer = APD -> Whitening filter ->
% matched filter -> equalizer (in continuous time)
H = Hw.*apd.H(sim.f).*rx.eq.Hrx.*rx.eq.Hff(sim.f/mpam.Rs).*exp(1j*2*pi*sim.f/mpam.Rs*grpdelay(rx.eq.num, rx.eq.den, 1));
h2 = real(ifft(ifftshift(H)));
h = h2(1:sim.Mct*floor(rx.eq.Ntaps/2));
hh = [h2(end-sim.Mct*floor(rx.eq.Ntaps/2):end); h];

hh = hh.*conj(hh); % |h(t)|^2
hh = hh/abs(sum(hh)); % normalize
 
BW = trapz(sim.f, abs(H).^2); % shot noise bandwidth
Ssh = BW*conv(hh, Ssh);
Ssh = circshift(Ssh, [-round(grpdelay(hh, 1, 1)) 0]); % remove delay due to equalizer

% Add thermal noise
Ssh = Ssh + rx.N0/2*trapz(sim.f, abs(Hw.*rx.eq.Hrx.*rx.eq.Hff(sim.f/mpam.Rs)).^2); % filter noise BW (includes noise enhancement penalty)

% Normalize and sample
Sshd = Ssh(floor(sim.Mct/2)+1:sim.Mct:end);
Sshd = Sshd(Ndisc+1:end-Ndisc)/(link_gain*Pmax)^2;

%% Calculate error probabilities using Gaussian approximation for each transmitted symbol
pe_gauss = zeros(mpam.M, 1); % symbol error probability for each level
dat = gray2bin(dataTX, 'pam', mpam.M); % fix index mapping
for k = 1:Nsymb
    sig = sqrt(Sshd(k));
    
    if dat(k) == mpam.M-1
        pe_gauss(dat(k)+1) = pe_gauss(dat(k)+1) + qfunc((yd(k)-Pthresh(end))/sig);
    elseif dat(k) == 0
        pe_gauss(dat(k)+1) = pe_gauss(dat(k)+1) + qfunc((Pthresh(1)-yd(k))/sig);
    else 
        pe_gauss(dat(k)+1) = pe_gauss(dat(k)+1) + qfunc((Pthresh(dat(k) + 1) - yd(k))/sig);
        pe_gauss(dat(k)+1) = pe_gauss(dat(k)+1) + qfunc((yd(k) - Pthresh(dat(k)))/sig);
    end
end

pe_gauss = real(pe_gauss)/Nsymb;

bergauss_levels = pe_gauss/log2(mpam.M); % this corresponds to p(error | given symbol)p(symbol)
bergauss = sum(pe_gauss)/log2(mpam.M);

%% AWGN Approximation:
% AWGN (including noise enhancement penalty)
% mpam = mpam.adjust_levels(tx.Ptx*link_gain, tx.rexdB);
hh0 = round(grpdelay(hh, 1, 1));
hhd_prev = hh(hh0:-sim.Mct:1);
hhd_post = hh(hh0:sim.Mct:end);
hhd = [hhd_prev(end:-1:2); hhd_post];
hh0 = length(hhd_prev);
hhd = hhd/abs(sum(hhd));
if isfield(rx, 'eq') && isfield(rx.eq, 'Ntaps')
    ber_awgn = mpam.ber_awgn(@(P) 1/(Pmax*link_gain)*noise_std(Pmax*link_gain*(P*hhd(hh0) + (sum(hhd)-hhd(hh0)))));
    % In calculating noise_std for AWGN channel, all other symbols in hh
    % are assumed to be the highest PAM level.
else
    ber_awgn = mpam.ber_awgn(@(P) 1/(Pmax*link_gain)*noise_std(Pmax*link_gain*P));
end