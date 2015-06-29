%% Calculate BER of amplified IM-DD system through montecarlo simulation
function ber = ber_soa_montecarlo(mpam, tx, fiber, soa, rx, sim)

% Frequency
f = sim.f/sim.fs;

% Rescale levels to desired power
rex = 10^(-abs(tx.rexdB)/10); % extinction ratio. Defined as Pmin/Pmaxl
link_gain = tx.kappa*soa.Gain*fiber.link_attenuation(tx.lamb)*rx.R;
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

% Random sequence
dataTX = randi([0 mpam.M-1], 1, sim.Nsymb);

% Modulated PAM signal in discrete-time
xd = Plevels(gray2bin(dataTX, 'pam', mpam.M) + 1);
xt = 1/tx.kappa*reshape(kron(xd, mpam.pshape(0:sim.Mct-1)).', sim.N, 1);

% Generate optical signal
[Et, ~] = optical_modulator(xt, tx, sim);

% Fiber propagation
Et = fiber.linear_propagation(Et, sim.f, tx.lamb);

% Amplifier
et = soa.amp(Et, sim.fs);

% Optical bandpass filter
eo = ifft(fft(et).*ifftshift(rx.optfilt.H(f)));

%% Direct detection and add thermal noise
%% Add shot noise
if isfield(sim, 'shot') && sim.shot
    q = 1.60217657e-19;      % electron charge (C)

    % Instataneous received power considering only attenuation from the fiber   
    Sshot = 2*q*(rx.R*abs(eo).^2 + rx.Id);     % one-sided shot noise PSD

    % Frequency is divided by two because PSD is one-sided
    wshot = sqrt(Sshot*sim.fs/2).*randn(size(eo));
else 
    wshot = 0;
end

% Direct detection and add noises
yt = abs(eo).^2;
yt = yt + wshot + sqrt(rx.N0*sim.fs/2)*randn(size(eo));

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