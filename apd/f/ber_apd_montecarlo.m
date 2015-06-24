function ber = ber_apd_montecarlo(mpam, tx, fiber, apd, rx, sim)

% time and frequency measures
f = sim.f/sim.fs;

% Random sequence
dataTX = randi([0 mpam.M-1], 1, sim.Nsymb);

% Rescale to desired power
rex = 10^(-abs(tx.rexdB)/10); % extinction ratio. Defined as Pmin/Pmax
link_gain = tx.kappa*fiber.link_attenuation(tx.lamb)*apd.R*apd.Gain;
if strcmp(mpam.level_spacing, 'uniform') 
    % Adjust levels to desired transmitted power and include additional dc bias due to finite extinction ratio
    Pmin = 2*tx.Ptx*rex/(1 + rex); % power of the lowest level 
    Plevels = mpam.a*(tx.Ptx/mean(mpam.a))*((1-rex)/(1+rex)) + Pmin; % levels at the transmitter
    
    Pthresh = mpam.b*(link_gain*tx.Ptx/mean(mpam.a))*((1-rex)/(1+rex)) + link_gain*Pmin; % decision thresholds at the #receiver#
elseif strcmp(mpam.level_spacing, 'nonuniform') % already includes extinction ratio penalty, so just scale
    % Adjust levels to desired transmitted power.
    % Extinction ratio penalty was already included in the level
    % optimization
    Plevels = mpam.a*tx.Ptx/mean(mpam.a); % levels at the transmitter
    
    Pthresh = mpam.b*link_gain*tx.Ptx/mean(mpam.a); % decision thresholds at the #receiver#
end  

% Modulated PAM signal in discrete-time
xd = Plevels(gray2bin(dataTX, 'pam', mpam.M) + 1);
xt = 1/tx.kappa*reshape(kron(xd, mpam.pshape(0:sim.Mct-1)).', sim.N, 1);

% Generate optical signal
[Et, ~] = optical_modulator(xt, tx, sim);

% Fiber propagation
[~, Pt] = fiber.linear_propagation(Et, sim.f, tx.lamb);

%% Detect and add noises
yt = apd.detect(Pt, sim.fs, 'gaussian', rx.N0);

% Electric low-pass filter
yt = real(ifft(fft(yt).*ifftshift(rx.elefilt.H(f))));

% Sample
yd = yt(sim.Mct/2:sim.Mct:end);

% Heuristic pdf for a level
if sim.verbose
    figure(100)
    [nn, xx] = hist(yd(dataTX == 2), 50);
    nn = nn/trapz(xx, nn);
    bar(xx, nn)
end

% Discard first and last sim.Ndiscard symbols
ndiscard = [1:sim.Ndiscard sim.Nsymb-sim.Ndiscard+1:sim.Nsymb];
yd(ndiscard) = []; 
dataTX(ndiscard) = [];

% Demodulate
dataRX = sum(bsxfun(@ge, yd, Pthresh.'), 2);
dataRX = bin2gray(dataRX, 'pam', mpam.M).';

% True BER
[~, ber] = biterr(dataRX, dataTX);
