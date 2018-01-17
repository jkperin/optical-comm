%% Validate calculations of GN model coefficients
clear, clc, close all

Fiber = fiber(50e3, @(l) 0.2);

lcenter = 1550e-9;
fcenter = lambda2freq(lcenter);
f = (-200:50:200)*1e9 + fcenter;
lamb = freq2lambda(f);
Df = 50e9;
Nspans = 100;

%% Test GN model
D = GN_model_coeff(lamb, Df, Fiber, 1);

figure, imagesc(10*log10(D))

