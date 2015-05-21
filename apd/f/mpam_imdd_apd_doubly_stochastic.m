function ber = mpam_imdd_apd_doubly_stochastic(mpam, tx, rx, sim)

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

% Transmitter pulse shape
gtx = ones(1, sim.Mct);

% Continuous time
xt = kron(x, gtx);

xt = reshape(xt.', 1, sim.Nsymb*sim.Mct);

% Scale so that average power is equal to the received power
xmean = mean(a); % this is the true mean
P = xt/xmean*tx.Prec;

% %% RIN noise
% nrin = 0;
% if sim.rin
%     % Calculate double-sided RIN PSD, which depends on the average power
%     Srin = 10^(tx.RIN/10)*(rx.Gapd*rx.R*P).^2;
% 
%     nrin = sqrt(Srin*2*mpam.bw).*randn(size(P));
% end

%% Shot noise
Ts = (log2(mpam.M)/mpam.Rb);
dt = Ts/sim.Mct;
df = 1/(length(P)*dt);
fs = 1/dt;
f = -fs/2:df:fs/2-df;

Irx = rx.Gapd*rx.R*P;
if true %strcmp(sim.shot, 'on')  
    Irx = apd_doubly_stochastic(P, dt, tx, rx); 
   
end

Idet = Irx + sqrt(rx.N0*mpam.Rs)*randn(size(P));

% Idet = filter(fliplr(gtx), 1, Idet);

% % Sample
% Idet = Idet(sim.Mct:sim.Mct:end);
% Idet = Irx(1:sim.Mct:end);
% Idet = Irx;

% Rescale signal for detection
y = xmean/mean(Idet)*Idet.';

% Demodulate
dataRX = sum(bsxfun(@ge, y, b.'), 2);
dataRX = bin2gray(dataRX, 'pam', mpam.M);

% True BER
Ndisc = 20;
[~, ber] = biterr(dataRX(Ndisc+1:end-Ndisc), dataTX(Ndisc+1:end-Ndisc));

% BER using approximation BER = SER/log2(M)
% ber = mean(dataTX ~= dataRX)/log2(mpam.M);
