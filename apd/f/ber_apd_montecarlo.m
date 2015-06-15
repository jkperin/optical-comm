function ber = ber_apd_montecarlo(mpam, tx, apd, rx, sim)

% time and frequency measures
f = sim.f/sim.fs;

% Random sequence
dataTX = randi([0 mpam.M-1], 1, sim.Nsymb);

% Modulated PAM signal in discrete-time
Pd = mpam.a(gray2bin(dataTX, 'pam', mpam.M) + 1);

% Add pulse 
Pt = reshape(kron(Pd, mpam.pshape(1:sim.Mct)).', sim.N, 1);

% Rescale to desired power
if strcmp(mpam.level_spacing, 'uniform') % adjust levels to include extinction ratio penalty
    Pmin = mean(Pt)/(10^(abs(tx.rex)/10)-1); % minimum power 
    Plevels = (mpam.a + Pmin)*apd.R*apd.Gain*tx.Ptx/(mean(Pt) + Pmin); % levels at the receiver
    Pthresh = (mpam.b + Pmin)*apd.R*apd.Gain*tx.Ptx/(mean(Pt) + Pmin); % decision thresholds at the receiver
    Pt = (Pt + Pmin)*tx.Ptx/(mean(Pt) + Pmin); % after rescaling E(Pt) = tx.Ptx
elseif strcmp(mpam.level_spacing, 'nonuniform') % already includes extinction ratio penalty, so just scale
    Plevels = mpam.a*apd.R*apd.Gain*tx.Ptx/mean(Pt);
    Pthresh = mpam.b*apd.R*apd.Gain*tx.Ptx/mean(Pt);
    Pt = Pt*tx.Ptx/mean(Pt);
end  

%% Add intensity noise
if isfield(sim, 'RIN') && sim.RIN
    % Calculate RIN PSD, which depends on the instantaneous power
    Srin = 10^(tx.RIN/10)*Pt.^2;

    % noise. Srin is two-sided psd
    wrin = sqrt(Srin.*sim.fs).*randn(size(Pt));

    Pt = Pt + wrin;
end 

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
Pd(ndiscard) = []; 
dataTX(ndiscard) = [];

% Demodulate
dataRX = sum(bsxfun(@ge, yd, Pthresh.'), 2);
dataRX = bin2gray(dataRX, 'pam', mpam.M).';

% True BER
[~, ber] = biterr(dataRX, dataTX);
