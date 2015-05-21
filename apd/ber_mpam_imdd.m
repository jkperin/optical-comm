function ber = ber_mpam_imdd(mpam, tx, rx, sim)

% Generate data
dataTX = randi([0 mpam.M-1], [sim.Nsymb 1]);

% Generate unipolar PAM signal
if strcmp(sim.levelSpacing, 'uniform')
    a = (0:2:2*(mpam.M-1)).';
    b = (1:2:(2*(mpam.M-1)-1)).';
elseif strcmp(sim.levelSpacing, 'nonuniform')
    a = mpam.a;
    b = mpam.b;
else
    error('Invalide Option!')
end

x = a(gray2bin(dataTX, 'pam', mpam.M) + 1);

% Scale so that average power is equal to the received power
xmean = mean(a); % this is the true mean
P = x/xmean*tx.Prec;

%% Thermal Noise
nthermal = sqrt(sim.N0*mpam.bw)*randn(size(P));

%% Detect and add noises
Idet = rx.detect(P, 1/mpam.bw) + nthermal;

% Rescale signal for detection
y = Idet/(rx.Gain*rx.R*tx.Prec)*xmean;

% Demodulate
dataRX = sum(bsxfun(@ge, y, b.'), 2);
dataRX = bin2gray(dataRX, 'pam', mpam.M);

% True BER
[~, ber] = biterr(dataRX, dataTX);

% BER using approximation BER = SER/log2(M)
% ber = mean(dataTX ~= dataRX)/log2(mpam.M);
