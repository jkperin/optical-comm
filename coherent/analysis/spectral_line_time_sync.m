%% Spectral line time recovery
clear, clc, close all

addpath ../../f/
addpath ../../mpam/

Nsymb = 16;
Mct = 61;
N = Nsymb*Mct;
Rb = 10e9;

pulse_shape = select_pulse_shape('rect', Mct);
mpam = PAM(2, Rb, 'equally-spaced', pulse_shape);
mpam = mpam.set_levels([-1; 1], 0);

fs = mpam.Rs*Mct;
[f, t] = freq_time(N, fs);

dataTX = randi([0 mpam.M-1], [1 Nsymb]);
dataTX(1:4) = [1 0 1 0];
x = mpam.signal(dataTX);
x = x - mean(mpam.a);

F = ClassFilter('bessel', 5, mpam.Rs/(fs/2), fs);
xf = F.filter(x);

BW = 1e9;
lpf = design_filter('butter', 5, BW/(fs/2));
[bpf.num, bpf.den] = iirlp2bp(lpf.num, lpf.den, BW/(fs/2), mpam.Rs/(fs/2) + BW/(fs/2)*[-1 1]); % converts to BPF
bpf.H = @(f) freqz(bpf.num, bpf.den, 2*pi*f);

x2 = xf.^2;

[d, w] = grpdelay(bpf.num, bpf.den);
figure, plot(fs/(2e9)*w/pi, d)
dd = interp1(fs/(2e9)*w/pi, d, mpam.Rs/1e9)

figure
plot(f/1e9, abs(fftshift(fft(x2))).^2)
plot(f/1e9, 1e6*abs(bpf.H(f/fs)).^2);
tdelay = 1/(4*mpam.Rs);
% xrec = real(ifft(fft(x2).*ifftshift(bpf.H(f/fs).*exp(-1j*2*pi*f*(tdelay - dd/fs)))));
xrec = filtfilt(bpf.num, bpf.den, x2);
figure, hold on
plot(xrec)
xrec = real(ifft(fft(xrec).*ifftshift(exp(-1j*2*pi*f*(tdelay - 0*dd/fs)))));
plot(xrec)

model = fit(t.', x2.'-mean(x2), 'sin1')
s = model(t').';

% Limiting amps
clk = xrec;
clk(clk > 0) = 1;
clk(clk < 0) = 0;

% Rising-edge detection
samp = diff([clk 0]);
idx = find(samp > 0);
clk = samp;
clk(clk < 0) = 0;

figure, hold on
K = 1:N;
plot(K, x)
plot(K, xf)
plot(K, x2)
% plot(t, s)
plot(K, 50*xrec, 'k')
plot(K, clk)


