%% Calculate BER of amplified IM-DD system through montecarlo simulation
function ber = ber_soa_montecarlo(mpam, tx, soa, rx, sim)

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
    Plevels = (mpam.a + Pmin)*soa.Gain*tx.Ptx/(mean(Pt) + Pmin); % levels at the receiver
    Pthresh = (mpam.b + Pmin)*soa.Gain*tx.Ptx/(mean(Pt) + Pmin); % decision thresholds at the receiver
    Pt = (Pt + Pmin)*tx.Ptx/(mean(Pt) + Pmin); % after rescaling E(Pt) = tx.Ptx
elseif strcmp(mpam.level_spacing, 'nonuniform') % already includes extinction ratio penalty, so just scale
    Plevels = mpam.a*soa.Gain*tx.Ptx/mean(Pt);
    Pthresh = mpam.b*soa.Gain*tx.Ptx/mean(Pt);
    Pt = Pt*tx.Ptx/mean(Pt);
end   

% Calculate electric field (no chirp) before amplifier
x = sqrt(Pt);

% Amplifier
et = soa.amp(x, sim.fs);
Ef = fftshift(fft(et));

% Optical bandpass filter
eo = ifft(fft(et).*ifftshift(rx.optfilt.H(f)));

% Direct detection and add thermal noise
yt = abs(eo).^2 + sqrt(rx.N0*sim.fs/2)*randn(size(eo));

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