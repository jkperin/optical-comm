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

% Automatic gain control
% Pmax = 2*tx.Ptx/(1 + 10^(-abs(tx.rexdB)/10)); % calculated from mpam.a
yt = yt/(Pmax*link_gain); % just refer power values back to transmitter
mpam = mpam.norm_levels;

%% Equalization
if isfield(rx, 'eq') && (isfield(tx, 'modulator') || ~isinf(apd.BW))
    rx.eq.type = strrep(rx.eq.type, 'Adaptive', 'Fixed'); % replace adaptive for fixed
    % Note: in this simulation there aren't many symbols to determine the
    % adaptive equalizer
else % otherwise only filter using rx.elefilt
    rx.eq.type = 'None';
end

% Equalizer
[yd, rx.eq] = equalize(rx.eq, yt, Hch, mpam, rx, sim);

% Symbols to be discard in BER calculation
yd = yd(Nzero+1:end-Nzero);

%% Detection
Pthresh = mpam.b; % refer decision thresholds to receiver

% Noise bandwidth
Df  = rx.elefilt.noisebw(sim.fs)/2;
% if strcmp(rx.eq.type, 'None')
%     Heq = 1;
% else
%     Heq = freqz(rx.eq.num, rx.eq.den, sim.f, mpam.Rs);
% end
% Htot = Hch.*Heq; 
% Df  = (1/(abs(Htot(sim.f==0)).^2))*trapz(sim.f, abs(Htot).^2)/2; % matched filter noise BW
% Dfshot = (1/(abs(apd.H(0).*Htot(sim.f==0)).^2))*trapz(sim.f, abs(apd.H(sim.f).*Htot).^2)/2;
% DfRIN = (1/(abs(fiber.H(0, tx).*apd.H(0).*Htot(sim.f==0)).^2))*trapz(sim.f, abs(fiber.H(sim.f, tx).*apd.H(sim.f).*Htot).^2)/2;

% Variance of thermal noise
varTherm = rx.N0*Df; 

if RIN
    varRIN =  @(Pnorm) 10^(tx.RIN/10)*(Pmax*link_gain*Pnorm).^2*Df;
else
    varRIN = @(Pnorm) 0;
end

%% Calculate error probabilities using Gaussian approximation for each transmitted symbol
pe_gauss = zeros(mpam.M, 1);
dat = gray2bin(dataTX, 'pam', mpam.M);
for k = 1:Nsymb
    mu = yd(k);
    varShot = apd.varShot(Pmax*link_gain*yd(k)/apd.Gain, Df); 
    % Note: apd.BW*pi/2 is the noise BW of the APD frequency response
%     sig = 1/(Pmax*link_gain)*sqrt(varTherm + varShot + varRIN(yd(k)));
    sig = 1/(Pmax*link_gain)*sqrt(rx.eq.Kne)*sqrt(varTherm + varShot + varRIN(yd(k)));
     
    if dat(k) == mpam.M-1
        pe_gauss(dat(k)+1) = pe_gauss(dat(k)+1) + qfunc((mu-Pthresh(end))/sig);
    elseif dat(k) == 0;
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