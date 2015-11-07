%% Calculate BER for unamplified IM-DD link with APD
% bergauss = BER using gaussian approximation
% bertail_levels = BER of each level individually. This actually
% corresponds bertail_levels = p(error | given symbol)p(symbol)

function [bergauss, bergauss_levels] = ber_apd_gauss(mpam, tx, fiber, apd, rx, sim)

Nsymb = mpam.M^sim.L; % number of data symbols 
% number of zero symbols pad to begin and end of sequnece.
if isfield(rx, 'eq') && isfield(rx.eq, 'Ntaps')
    Nzero = max(rx.eq.Ntaps, sim.L);
else
    Nzero = sim.L;
end
N = sim.Mct*(Nsymb + 2*Nzero); % total number of points 

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
xt = mpam.mod(dataTX, sim.Mct);
xt = [zeros(sim.Mct*Nzero, 1); xt; zeros(sim.Mct*Nzero, 1)]; % zero pad

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

%% Automatic gain control
% Normalize signal so that highest level is equal to 1
yt = yt/(Pmax*link_gain); 
mpam = mpam.norm_levels;

%% Equalization
[yd, rx.eq] = equalize(rx.eq, yt, Hch, mpam, rx, sim);

% Symbols to be discard in BER calculation
yd = yd(Nzero+1:end-Nzero);

%% Detection
Pthresh = mpam.b; % refer decision thresholds to receiver

% Noise bandwidth
noise_std = apd.stdNoise(rx.eq.Hrx, rx.eq.Hff(sim.f/mpam.Rs), rx.N0, tx.RIN, sim);
    
%% Calculate error probabilities using Gaussian approximation for each transmitted symbol
pe_gauss = zeros(mpam.M, 1);
dat = gray2bin(dataTX, 'pam', mpam.M);
for k = 1:Nsymb
    mu = yd(k);
    sig = 1/(Pmax*link_gain)*noise_std(Pmax*link_gain*yd(k));
     
    if dat(k) == mpam.M-1
        pe_gauss(dat(k)+1) = pe_gauss(dat(k)+1) + qfunc((mu-Pthresh(end))/sig);
    elseif dat(k) == 0
        pe_gauss(dat(k)+1) = pe_gauss(dat(k)+1) + qfunc((Pthresh(1)-mu)/sig);
    else 
        pe_gauss(dat(k)+1) = pe_gauss(dat(k)+1) + qfunc((Pthresh(dat(k) + 1) - mu)/sig);
        pe_gauss(dat(k)+1) = pe_gauss(dat(k)+1) + qfunc((mu - Pthresh(dat(k)))/sig);
    end
end

pe_gauss = real(pe_gauss)/Nsymb;

bergauss_levels = pe_gauss/log2(mpam.M); % this corresponds to p(error | given symbol)p(symbol)
bergauss = sum(pe_gauss)/log2(mpam.M);

if sim.verbose
    figure(102), hold on
    plot(link_gain*xt)
    plot(ix, yd, 'o')
    legend('Transmitted power', 'Sampled received signal')
end