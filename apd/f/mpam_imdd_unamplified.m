function ber = mpam_imdd_unamplified(mpam, tx, rx, sim)

% Generate data
dataTX = randi([0 mpam.M-1], [sim.Nsymb 1]);

% Generate unipolar PAM signal
if strcmp(sim.levelSpacing, 'uniform')
    a = (0:2:2*(mpam.M-1)).';
    b = (1:2:(2*(mpam.M-1)-1)).';
elseif strcmp(sim.levelSpacing, 'nonuniform')
    [a, b] = calcOptLevelSpacing(mpam, tx, rx, sim);   
else
    error('Invalide Option!')
end

x = a(gray2bin(dataTX, 'pam', mpam.M) + 1);

a/a(2);

% Scale so that average power is equal to the received power
xmean = mean(a); % this is the true mean
P = x/xmean*tx.Prec;

%% RIN noise
nrin = 0;
if strcmp(sim.rin, 'on')
    % Calculate double-sided RIN PSD, which depends on the average power
    Srin = 10^(tx.RIN/10)*(rx.Gapd*rx.R*P).^2;

    nrin = sqrt(Srin*2*mpam.bw).*randn(size(P));
end

%% Shot noise
nshot = 0;
if strcmp(sim.shot, 'on')
    q = 1.60217657e-19;      % electron charge (C)
    Id = 0;                  % dark current

    % Instataneous received power considering only attenuation from the fiber   
    Sshot = rx.Gapd^2*rx.R^2*rx.Fa(rx.ka, rx.Gapd)*2*q*(rx.R*P + Id);     % one-sided shot noise PSD

    % Frequency is divided by two because PSD is one-sided
    nshot = sqrt(Sshot*mpam.bw).*randn(size(P));
end

%% Thermal Noise
nthermal = sqrt(rx.N0*mpam.bw)*randn(size(P));

%% Detect and add noises
Idet = rx.Gapd*rx.R*P + nrin + nshot + nthermal;

% Rescale signal for detection
y = Idet/(rx.Gapd*rx.R*tx.Prec)*xmean;

% Demodulate
dataRX = sum(bsxfun(@ge, y, b.'), 2);
dataRX = bin2gray(dataRX, 'pam', mpam.M);

% True BER
[~, ber] = biterr(dataRX, dataTX);

% BER using approximation BER = SER/log2(M)
% ber = mean(dataTX ~= dataRX)/log2(mpam.M);
