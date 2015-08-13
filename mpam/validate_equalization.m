%% Compare performance of M-PAM in 3 different scenarios
% 1. PIN receiver, no amplifier
% 2. SOA & PIN receiver
% 3. APD
clear, clc, close all

addpath ../f/

%% Simulation parameters
sim.Nsymb = 2^7; % Number of symbols in montecarlo simulation
sim.Mct = 45;     % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
sim.BERtarget = 1e-4; 
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

sim.ros = 2;

%% M-PAM
mpam = PAM(4, 100e9, 'equally-spaced', @(n) double(n >= 0 & n < sim.Mct));

%% Time and frequency
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'

dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

sim.t = t;
sim.f = f;

% Modulator frequency response
tx.modulator.fc = 50e9; % modulator cut off frequency
tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds
tx.modulator.H = @(f) (1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2)).*exp(1j*2*pi*f*tx.modulator.grpdelay);
tx.modulator.h = @(t) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0));

%% Receiver
rx.N0 = (20e-12).^2; % thermal noise psd
rx.Id = 10e-9; % dark current
rx.R = 1; % responsivity
% Electric Lowpass Filter
% rx.elefilt = design_filter('bessel', 5, mpam.Rs/(sim.fs/2));
rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);
rx.antialiasing = design_filter('bessel', 4, mpam.Rs/(sim.fs/2));

dataTX = randi([0 mpam.M-1], 1, sim.Nsymb);
x = mpam.mod(dataTX, sim.Mct);

p = real(ifft(fft(x).*ifftshift(tx.modulator.H(f))));
% p =x;

% p = p + sqrt(rx.N0*sim.fs/2)*randn(size(p));

pa = real(ifft(fft(p).*ifftshift(rx.antialiasing.H(f/sim.fs))));

figure, hold on
plot(t, x)
plot(t, pa)

%% Time-domain fractionally-spaced equalizer

% pk = pa(1:floor(sim.Mct/sim.ros):end);
% tk = t(1:floor(sim.Mct/sim.ros):end);

pk = resample(pa, sim.ros, sim.Mct);
tk = resample(sim.t, sim.ros, sim.Mct);

plot(tk, pk, 'o')

Ntaps = 7;
Ntrain = 10e3;
W = [zeros(Ntaps-1, 1); 1];
y = zeros(size(pk));
b = mpam.mod(dataTX, 1);
mu = 1e-2;
n = 1;
for k = Ntaps:length(pk)
    z = pk(k-Ntaps+1:k);
    
    y(k) = sum(W.*z);
    
    if mod(k, sim.ros) == 0
        if n < Ntrain % Training
%             e((k-1)/sim.ros) = y(k) - b((k-1)/sim.ros);

            e(k/sim.ros) = y(k) - b(k/sim.ros-1);
        
        else
%             e((k-1)/sim.ros) = y(k) - mpam.mod(mpam.demod(y(k)), 1);
            
            e(k/sim.ros) = y(k) - mpam.mod(mpam.demod(y(k)), 1);
            
        end
        n = n + 1;
        
%         W = W - 2*mu*e((k-1)/sim.ros)*z; 
        W = W - 2*mu*e(k/sim.ros)*z;
    end
end

% plot(tk, b, 'o')

plot(tk(2:sim.ros:end), y(2:sim.ros:end), 'o')
figure, plot(e)
        
[H, w] = freqz(W(end:-1:1), 1);
figure, hold on
plot(2*mpam.Rs/1e9*w/(2*pi), abs(H).^2)
% ff = 2*mpam.Rs*w/(2*pi);
% X = tx.modulator.H(sim.f).*rx.antialiasing.H(sim.f/sim.fs).*rx.elefilt.H(sim.f/sim.fs);
% x = real(ifft(ifftshift(X)));
% x = x(1:sim.Mct:end);
% Xf = fftshift(fft(x))*sim.Mct;
% Heq = 1./(Xf);
% plot(abs(Heq).^2)
        
%% Time-domain symbol-rate equalizer
% matchedfilt = design_filter('matched', @(t) conv(mpam.pshape(t), 1/sim.fs*tx.modulator.h(t/sim.fs), 'full') , 1/sim.Mct); 
matchedfilt = design_filter('matched', mpam.pshape, 1/sim.Mct); 

Hmatched = conj(tx.modulator.H(sim.f).*matchedfilt.H(sim.f/sim.fs));

% ps = real(ifft(fft(p).*ifftshift(matchedfilt.H(f/sim.fs))));
ps = real(ifft(fft(p).*ifftshift(Hmatched)));

ps = ps(floor(sim.Mct/2):sim.Mct:end);
ts = t(floor(sim.Mct/2):sim.Mct:end);

Ntaps = 15;
Ntrain = 10e3;
Ws = [zeros(Ntaps-1, 1); 1];
ys = zeros(size(ps));
mu = 1e-2;
n = 1;
for k = Ntaps:length(ps)
    z = ps(k-Ntaps+1:k);
    
    ys(k) = sum(Ws.*z);
    
    if n < Ntrain % Training
%             e((k-1)/sim.ros) = y(k) - b((k-1)/sim.ros);

        es(k) = ys(k) - b(k);

    else
%             e((k-1)/sim.ros) = y(k) - mpam.mod(mpam.demod(y(k)), 1);

        es(k) = ys(k) - mpam.mod(mpam.demod(ys(k)), 1);

    end
    n = n + 1;

%         W = W - 2*mu*e((k-1)/sim.ros)*z; 
    Ws = Ws - 2*mu*es(k)*z;
end

% figure, plot(ts, ys, 'o')
figure, plot(es)

[Hs, fss] = freqz(Ws(end:-1:1), 1, 200, mpam.Rs);
figure, hold on
plot(fss/1e9, abs(Hs).^2)

X = abs(Hmatched).^2;
X = conj(flipud(X(1:length(X)/2+1)));
ff = -flipud(sim.f(1:length(sim.f)/2+1));
Nfold = floor(mpam.Rs/(2*df))+1;
Xfolded = X(1:Nfold);
ffolded = ff(1:Nfold);
fold = true;
for k = Nfold+1:Nfold:length(ff)-Nfold
    if fold
        [ff(k+Nfold-1) ff(k)]/1e9
        
        Xfolded = Xfolded + conj(flipud(X(k:k+Nfold-1)));
        
        fold = false;
    else
        [ff(k) ff(k+Nfold-1)]/1e9
        
        Xfolded = Xfolded + X(k:k+Nfold-1);

        fold = true;
    end
end

Xfolded = [flipud(conj(Xfolded(2:end))); Xfolded(1:end-1)];
ffolded = [-flipud((ffolded(2:end))); ffolded(1:end-1)];

x = ifft(ifftshift(abs(matchedfilt.H(sim.f/sim.fs).*tx.modulator.H(sim.f)).^2));
x = x(1:sim.Mct:end);
Xfolded2 = fftshift(fft(x))*sim.Mct;

plot(ffolded/1e9, 1./Xfolded)
plot(ffolded/1e9, 1./Xfolded2)
legend('Adaptive', 'Theory folded', 'Theory downsampled')

pf = real(ifft(fft(ps)./ifftshift(Xfolded)));

figure, hold on
plot(ts, ps, '*')
plot(ts, ys, 's')
plot(ts, pf, 'o')
legend('sampled', 'Adaptive', 'Theory')




%% Equalizer
% Nfft = sim.Nsymb*sim.ros;
% ff = -0.5:1/Nfft:0.5-1/Nfft;
% ef = design_filter('matched', mpam.pshape, 1/sim.ros);
% % W = 1./(tx.modulator.H(ff*mpam.Rs*sim.ros).*ef.H(ff));
% W = 1;
% % figure
% % plot(
% 
% Pk = fftshift(fft(pk));
% Pk = W.'.*Pk;
% 
% Ns = round((sim.ros-1)/sim.ros*length(pk))/2;
% % Pk = [0:length(pk)/2 length(pk)/2-1:-1:1];
% 
% Pk2 = Pk(Ns+1:end-Ns);
% Pk2(1:Ns+1) = Pk2(1:Ns+1) + conj(foldlr(Pk(1:Ns+1)));
% Pk2(end-Ns+2:end) = Pk2(end-Ns+2:end) + conj(foldlr(Pk(end-Ns+2:end)));
% 
% figure, hold on
% plot(-128:127, abs(Pk))
% plot(-64:63, abs(Pk2))
% 
% yk = real(ifft(ifftshift(Pk2)));
% 
% figure(1)
% plot(tk(1:sim.ros:end), yk/2, '-*')


