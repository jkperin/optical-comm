%% Analog duobinary filter candidate impulse response
clear, clc, close all

addpath ../f/

Rs = 10e9;
fs = 200e9;

filt = design_filter('Bessel', 5, 1.8e9/(fs/2));

t = 0:100;
h = impz(filt.num, filt.den, t);
delay = 0*grpdelay(filt.num, filt.den, 1);
h = interp1(t-delay, h, t);
tn = -fs/Rs*5:fs/Rs:fs/Rs*5;
hn = interp1(t, h, tn);

figure, hold on
plot(t, h)
stem(tn, hn)
