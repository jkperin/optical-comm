%% higher-order Gaussian filters have a large negative excursion

clear, clc, close all

%% Gaussian vs Bessel

B_Gaussian = 2;
N_Gaussian = 4;

B_Bessel = 2;
N_Bessel = 5;

B_Butter = 2;
N_Butter = 5;


N = 1024;

fs = 20;
df = fs/N;
f = -fs/2:df:fs/2-df;

dt = 1/(fs);
t = -N/2*dt:dt:(N/2-1)*dt;

% Gaussian
F_0 = B_Gaussian/(2*log(2)^(1/(2*N_Gaussian)));

H_gaussian = exp(-0.5*(f/F_0).^(2*N_Gaussian));
H_gaussian = H_gaussian/max(abs(H_gaussian));

h_gaussian = fftshift(ifft(ifftshift(H_gaussian)));

% Bessel
F_Cut = B_Bessel/2;                   %% Normalise cut-off frequency

W_0 = 2*pi*F_Cut/0.615685;                               %% Constant group delay frequency for "besself" function (heuristic)

[b,a] = besself(N_Bessel, W_0);                             %% Design filter using Matlab tool (Poles and Zeros in S-domain)

H_bessel = freqs(b, a, 2*pi*f);
H_bessel = H_bessel/max(abs(H_bessel));

h_bessel = fftshift(ifft(ifftshift(H_bessel)));

% Butter
F_cut = B_Butter/2;
[b,a] = butter(N_Butter, 2*pi*F_cut, 's');
H_Butter = freqs(b, a, 2*pi*f);
H_Butter = H_Butter/max(abs(H_Butter));
h_Butter = fftshift(ifft(ifftshift(H_Butter)));


%% Figure
figure
plot(f, abs(H_gaussian).^2, 'k', f, abs(H_bessel).^2, 'r', f, abs(H_Butter), 'b')
xlabel('Frequency')
ylabel('|H(f)|^2')

figure
plot(t, h_gaussian, 'k', t, h_bessel, 'r', t, h_Butter, 'b')
axis([-5 5 -0.1 0.5])


Mct = 8;
cutoff = 1/1.23;
order = 1;

nlength = 256;   % makes it even, so grp delay is (tx.interp_length-1)/2

% Gaussian filter
f = linspace(-0.5, 0.5, nlength);
f0 = (cutoff/Mct)/(2*log(2)^(1/(2*order)));
H = exp(-0.5*(f/f0).^(2*order));

h = fftshift(ifft(ifftshift(H)));

afilt = 1;
bfilt = h;

H2 = freqz(bfilt, afilt, f, 1);

figure
plot(f, abs(H),f, abs(H2))
