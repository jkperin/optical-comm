function ber = ber_apd_montecarlo(mpam, tx, fiber, apd, rx, sim)
verbose = sim.verbose;
sim.verbose = max(sim.verbose-1, 0);

%% Channel response
% Hch does not include transmitter or receiver filter
if isfield(tx, 'modulator')
    Hch = tx.modulator.H(sim.f).*exp(1j*2*pi*sim.f*tx.modulator.grpdelay)...
    .*fiber.H(sim.f, tx).*apd.H(sim.f);
else
    Hch = fiber.H(sim.f, tx).*apd.H(sim.f);
end

link_gain = apd.Gain*apd.R*fiber.link_attenuation(tx.lamb); % Overall link gain

% Ajust levels to desired transmitted power and extinction ratio
mpam = mpam.adjust_levels(tx.Ptx, tx.rexdB);
Pmax = mpam.a(end); % used in the automatic gain control stage

% Modulated PAM signal
dataTX = randi([0 mpam.M-1], 1, sim.Nsymb); % Random sequence
xt = mpam.mod(dataTX, sim.Mct);
xt(1:sim.Mct*sim.Ndiscard) = 0; % zero sim.Ndiscard first symbols
xt(end-sim.Mct*sim.Ndiscard+1:end) = 0; % zero sim.Ndiscard last symbbols

% Generate optical signal
[Et, ~] = optical_modulator(xt, tx, sim);

% Fiber propagation
[~, Pt] = fiber.linear_propagation(Et, sim.f, tx.lamb);

%% Detect and add noises
yt = apd.detect(Pt, sim.fs, 'gaussian', rx.N0);

%% Whitening filter
if sim.WhiteningFilter
    [Hw, yt] = apd.Hwhitening(sim.f, tx.Ptx, rx.N0, yt);
else
    Hw = 1;
end

%% Automatic gain control
% Normalize signal so that highest level is equal to 1
yt = yt/(Pmax*link_gain);
% shot = shot/(Pmax*link_gain);
mpam = mpam.norm_levels;

%% Equalization
[yd, rx.eq] = equalize(rx.eq, yt, Hw.*Hch, mpam, rx, sim);

% Symbols to be discarded in BER calculation
Ndiscard = sim.Ndiscard*[1 1];
if isfield(rx.eq, 'Ntrain')
    Ndiscard(1) = Ndiscard(1) + rx.eq.Ntrain;
end
if isfield(rx.eq, 'Ntaps')
    Ndiscard = Ndiscard + rx.eq.Ntaps;
end
ndiscard = [1:Ndiscard(1) sim.Nsymb-Ndiscard(2):sim.Nsymb];
yd(ndiscard) = []; 
dataTX(ndiscard) = [];

% Demodulate
dataRX = mpam.demod(yd);

% True BER
[~, ber] = biterr(dataRX, dataTX);

if verbose  
    % Heuristic pdf for a level
    figure(100)
    [nn, xx] = hist(yd(dataTX == 2), 50);
    nn = nn/trapz(xx, nn);
    bar(xx, nn)
end    
