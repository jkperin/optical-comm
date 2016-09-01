%% Generate Duobinary eye diagram
clear, clc, close all

addpath ../../f/
addpath ../../mpam/

Mct = 15;
Nsymb = 2^12;
N = Nsymb*Mct;
pulse_shape = select_pulse_shape('rect', Mct);
mpam = PAM(4, 56e9, 'equally-spaced', pulse_shape);
fs = mpam.Rs*Mct;
f = freq_time(N, fs);

mpamdb = mpam.set_levels(0:mpam.M-1, 0.5 + (0:mpam.M-2));    

dataTX = randi([0 mpam.M-1], [1 Nsymb]);
xd = mpamdb.mod(dataTX); % Modulated PAM signal

xd = duobinary_encoding(xd);

xd_enc = xd;

ximp = upsample(xd_enc, mpam.pulse_shape.sps);
xt = filter(ones(1, mpam.pulse_shape.sps), 1, ximp);

filt = design_filter('bessel', 5, mpam.Rs/(fs/2));

xt = xt  + 0.1*randn(size(xt));
xt = real(ifft(fft(xt).*ifftshift(filt.H(f/fs).*exp(1j*2*pi*f/fs*Mct/2))));


figure, hold on, box on
Ntraces = 300;
x = xt(100*Mct+1:Ntraces*Mct);
n = 2*Mct;
M = floor(length(x)/n);
sel = 1:M*n;

x = reshape(x(sel), n, M);
plot(1:n, x, 'k', 'LineWidth', 0.5)

a = axis;
axis([1 n a(3) a(4)])


m = matlab2tikz(gca);
m.write('duobinary-4pam-eyediagram.tex')


xxt = mpam.signal(dataTX);
xxt = xxt - mean(xxt);
xxt = real(ifft(fft(xxt).*ifftshift(filt.H(f/fs).*exp(1j*2*pi*f/fs*Mct/2))));

xt = xt/sqrt(mean(abs(xt).^2));
xxt = xxt/sqrt(mean(abs(xxt).^2));

figure, hold on, box on
% plot(f/1e9, abs(fftshift(fft(xt))).^2, 'k')
plot(f/1e9, 5e6*abs(sinc(f/28e9)).^2, 'r', 'LineWidth', 4)
plot(f/1e9, 5e6*abs(sinc(f/14e9)).^2, 'k', 'LineWidth', 4)
xlabel('Frequency (GHz)', 'FontSize', 18)
ylabel('Power spectrum', 'FontSize', 18)
set(gca, 'FontSize', 18)
set(gca, 'ytick', [])
axis([-30 30 0 5e6])
legend('Ordinary 4-PAM', 'Duobinary 4-PAM')

figure, hold on, box on
plot(f/1e9, abs(fftshift(fft(xxt))).^2, 'k')
plot(f/1e9, 5e6*abs(sinc(f/28e9)).^2, 'k', 'LineWidth', 4)
xlabel('Frequency (GHz)', 'FontSize', 18)
ylabel('Power spectrum', 'FontSize', 18)
set(gca, 'FontSize', 18)
set(gca, 'ytick', [])
axis([-30 30 0 8e6])