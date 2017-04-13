%% Number of taps for given fiber length
clear, clc, close all

addpath ../f

Lkm = 1:20;
f = fiber(0);
lamb = 1250e-9;
Rs = 56e9;
ros = 5/4;

for k = 1:length(Lkm)
    f.L = 1e3*Lkm(k);
    N(k) = ceil(f.Ntaps(Rs, ros, lamb));
    D(k) = f.D(lamb)*f.L;
end

figure, plot(Lkm, N)
figure, plot(abs(D), N)