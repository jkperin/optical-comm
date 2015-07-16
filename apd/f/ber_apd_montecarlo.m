function ber = ber_apd_montecarlo(mpam, tx, fiber, apd, rx, sim)

% Normalized frequency
f = sim.f/sim.fs;

% Overall link gain
link_gain = tx.kappa*apd.Gain*fiber.link_attenuation(tx.lamb)*rx.R;

% Ajust levels to desired transmitted power and extinction ratio
mpam.adjust_levels(tx.Ptx, tx.rexdB);

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

% Electric low-pass filter
yt = real(ifft(fft(yt).*ifftshift(rx.elefilt.H(f))));

% Sample
ix = (sim.Mct-1)/2+1:sim.Mct:length(yt); % sampling points
yd = yt(ix);

% Discard first and last sim.Ndiscard symbols
ndiscard = [1:sim.Ndiscard sim.Nsymb-sim.Ndiscard+1:sim.Nsymb];
yd(ndiscard) = []; 
dataTX(ndiscard) = [];

% Automatic gain control
yd = yd/link_gain; % just refer power values back to transmitter

% Demodulate
dataRX = mpam.demod(yd);

% True BER
[~, ber] = biterr(dataRX, dataTX);

if sim.verbose
    % Signal
    figure(102), hold on
    plot(link_gain*Pt)
    plot(yt, '-k')
    plot(ix, yd, 'o')
    legend('Transmitted power', 'Received signal', 'Samples')
    
    % Heuristic pdf for a level
    figure(100)
    [nn, xx] = hist(yd(dataTX == 2), 50);
    nn = nn/trapz(xx, nn);
    bar(xx, nn)
end    
