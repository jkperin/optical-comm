%% Calculate BER of amplified IM-DD system through montecarlo simulation
function ber = ber_soa_montecarlo(mpam, tx, fiber, soa, rx, sim)

% Normalized frequency
f = sim.f/sim.fs;

% Overall link gain
link_gain = soa.Gain*fiber.link_attenuation(tx.lamb)*rx.R;

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
Et = fiber.linear_propagation(Et, sim.f, tx.lamb);

% Amplifier
et = soa.amp(Et, sim.fs);

% Optical bandpass filter
eo = [ifft(fft(et(:, 1)).*ifftshift(rx.optfilt.H(f))),...
    ifft(fft(et(:, 2)).*ifftshift(rx.optfilt.H(f)))];

%% Direct detection and add thermal noise
%% Shot noise
if isfield(sim, 'shot') && sim.shot
    q = 1.60217657e-19;      % electron charge (C)

    % Instataneous received power considering only attenuation from the fiber   
    Sshot = 2*q*(rx.R*sum(abs(eo).^2,2)  + rx.Id);     % one-sided shot noise PSD

    % Frequency is divided by two because PSD is one-sided
    wshot = sqrt(Sshot*sim.fs/2).*randn(size(Et));
else 
    wshot = 0;
end

% Direct detection and add noises
if isfield(sim, 'polarizer') && ~sim.polarizer 
    yt = abs(eo(:, 1)).^2 + abs(eo(:, 2)).^2;
else % by default assumes that polarizer is used
    yt = abs(eo(:, 1)).^2;
end
yt = yt + wshot + sqrt(rx.N0*sim.fs/2)*randn(size(Et));

% Automatic gain control
% Pmax = 2*tx.Ptx/(1 + 10^(-abs(tx.rexdB)/10)); % calculated from mpam.a
yt = yt/(Pmax*link_gain); % just refer power values back to transmitter
mpam.norm_levels;

%% Equalization
if isfield(rx, 'eq')
    rx.eq.TrainSeq = dataTX;
else % otherwise only filter using rx.elefilt
    rx.eq.type = 'None';
end
   
% Equalize
[yd, rx.eq] = equalize(rx.eq.type, yt, mpam, tx, fiber, rx, sim);

% Symbols to be discard in BER calculation
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
