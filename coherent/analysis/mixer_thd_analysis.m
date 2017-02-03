%% Mixer THD analysis
% Calculates total harmonic distortion at the mixer output
clear, clc, close all

addpath ../analog/
addpath ../../f/

N = 10000;
fs = 100e9;
dt = 1/fs;
t = 0:dt:(N-1)*dt;
f0 = 0.5e9;
df = 1/(dt*N);
f = -fs/2:df:fs/2-df;

filt = design_filter('bessel', 3, 0.7*28e9/(fs/2));
N0 = 0;

M = AnalogMixer(filt, N0, fs);

x = sin(2*pi*f0*t);

y = M.mix(x, x);

X = fftshift(fft(x))/N;
Y = fftshift(fft(y))/N;

figure, hold on
% plot(f, abs(X))
plot(f, abs(Y))

p0 = find(f == 0);
pf0 = find(f == f0);
p2f0 = find(f == 2*f0);
p4f0 = find(f == 4*f0);

THD = 10*log10(sum(abs(Y(p4f0:100:end)).^2)/abs(Y(p2f0).^2))
THDm = thd(y)

THDp = 10^(THDm/10)

Vamp = 1:0.1:2; % Controls amount of distortion in the mixer

for k = 1:length(Vamp)
    M.Vamp = Vamp(k);
    M = M.reset();
    y = M.mix(x, x);
    
    THDl(k) = thd(y);
end

figure
plot(Vamp, THDl, '-ko')
title('Mixer THD @ 1 GHz')
xlabel('Vamp')
ylabel('THD (dBc)')

figure
plot(Vamp, 100*10.^(THDl/10), '-ko')
title('Mixer THD @ 1 GHz')
xlabel('Vamp')
ylabel('THD (%)')

