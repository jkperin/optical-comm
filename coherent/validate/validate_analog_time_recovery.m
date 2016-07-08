%% Validate analog_time_recovery.m
clear, clc, close all

addpath ../analog
addpath ../../f/
addpath ../../mpam/

sim.Nsymb = 2^12;
sim.Mct = 10;
sim.N = sim.Mct*sim.Nsymb;

%% M-PAM
pulse_shape = select_pulse_shape('rrc', sim.Mct, 0.25, 6);

mpam = PAM(2, 56e9, 'equally-spaced', pulse_shape);

% Time and frequency
sim.fs = mpam.Rs*sim.Mct;
df = sim.fs/sim.N;
f = -sim.fs/2:df:sim.fs/2-df;

sim.f = f;
sim.Rs = mpam.Rs;

%% Time recovery
% Squarer
TimeRec.type = 'squarer';
BW = 1e9;
lpf = design_filter('bessel', 5, BW/(sim.fs/2));
[bpf.num, bpf.den] = iirlp2bp(lpf.num, lpf.den, BW/(sim.fs/2), mpam.Rs/(sim.fs/2) + BW/(sim.fs/2)*[-1 1]);
TimeRec.bpf = bpf;
TimeRec.bpf.H = @(f) freqz(bpf.num, bpf.den, 2*pi*f).*exp(1j*2*pi*f*grpdelay(bpf.num, bpf.den, 1));

figure, hold on
Hlpf = freqz(lpf.num, lpf.den, f, sim.fs);
plot(f/1e9, abs(TimeRec.bpf.H(sim.f/sim.fs)).^2, f/1e9, abs(Hlpf).^2, '--')
a = axis;
axis([-100 100 a(3:4)]);
legend('BPF', 'LPF')

% M&M
TimeRec.type = 'Mueller-Muller';
TimeRec.csi = sqrt(2)/2;
TimeRec.wn = 2*pi*3e9;
TimeRec.CT2DT = 'bilinear';
TimeRec.detect = @(x) sign(x);

%% Generate signal, add noise and delay
dataTX = randi([0 1], [1 sim.Nsymb]);

x = mpam.signal(dataTX);

x = 2*(x - mean(x));

% add noise
x = x + 0.1*randn(size(x));

% filter
lpf = design_filter('butter', 5, 0.5*mpam.Rs/(sim.fs/2));
x = filter(lpf.num, lpf.den, x);

% delay 
x = circshift(x, [0 -20]);

% Time recovery
[ysamp, idx] = analog_time_recovery(x, TimeRec, sim);

% Plots
figure, eyediagram(x, 2*sim.Mct)
title('Received signal eye diagram')

figure, 
subplot(211)
plot(ysamp, 'o')
title(sprintf('Sampling done with %s time recovery method', TimeRec.type))
subplot(212)
plot(x((sim.Mct)/2:sim.Mct:end), 'o')
title('Blind sampling')

figure,
plot(diff(idx), 'o')

figure
ind = -20:20;
c = xcorr(mpam.mod(dataTX), ysamp, 20, 'coeff');
plot(ind, abs(c))

[~, p] = max(abs(c));

ysamp = circshift(ysamp, [0 ind(p)]);
mpam.b = 0;
biterr(mpam.demod(ysamp(100:end-100)), dataTX(100:end-100))



