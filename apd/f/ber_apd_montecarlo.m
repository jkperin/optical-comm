function ber = ber_apd_montecarlo(mpam, tx, fiber, apd, rx, sim)

% Normalized frequency
f = sim.f/sim.fs;

% Overall link gain
link_gain = apd.Gain*fiber.link_attenuation(tx.lamb)*apd.R;

% Ajust levels to desired transmitted power and extinction ratio
mpam.adjust_levels(tx.Ptx, tx.rexdB);
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

% Automatic gain control
% Pmax = 2*tx.Ptx/(1 + 10^(-abs(tx.rexdB)/10)); % calculated from mpam.a
yt = yt/(Pmax*link_gain); % just refer power values back to transmitter
mpam.norm_levels;

%% Equalization
if isfield(rx, 'eq')
    rx.eq.TrainSeq = dataTX;
    [yd, Heq] = equalize(rx.eq.type, yt, mpam, tx, fiber, rx, sim);
else % otherwise only filter using rx.elefilt
    yd = equalize('None', yt, mpam, tx, fiber, rx, [], sim);
end
   
% Discard first and last sim.Ndiscard symbols
ndiscard = [1:sim.Ndiscard sim.Nsymb-sim.Ndiscard+1:sim.Nsymb];
yd(ndiscard) = []; 
dataTX(ndiscard) = [];

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
