%% Validate calculations of GN model coefficients
clear, close all

Fiber = fiber(50e3, @(l) 0.2);

lcenter = 1550e-9;
fcenter = lambda2freq(lcenter);
f = (-200:50:200)*1e9 + fcenter;
lamb = freq2lambda(f);
Df = 50e9;
Nspans = 100;

%% Test GN model
% D = GN_model_coeff(lamb, Df, Fiber, 0);
% 
% figure, imagesc(10*log10(D))

%% Load from file
df = 50e9;
dlamb = df2dlamb(df);
lamb = 1522e-9:dlamb:1582e-9;

S = load('GN_model_coeff_spanLengthkm=50.mat');
nonlinear_coeff = S.nonlinear_coeff;

for n = 1:length(nonlinear_coeff)
    figure, imagesc(10*log10(nonlinear_coeff{n}))
    axis square
    colorbar    
end

P = dBm2Watt(0.25 + 2*(rand(size(lamb))-0.5));

tic
NL = GN_model_noise(P, nonlinear_coeff);
toc, tic
NL_m = GN_model_noise_m(P, nonlinear_coeff);
toc

% Scale NL noise to "Namp" spans
epsilon = 0.05; % From Fig. 17 of P. Poggiolini and I. Paper, “The GN Model
% of Non-Linear Propagation in Uncompensated Coherent Optical Systems,” 
% J. Lightw. Technol., vol. 30, no. 24, pp. 3857–3879, 2012.
NL = NL*40^(1+epsilon);
NL_m = NL_m*40^(1+epsilon);


figure, hold on
plot(lamb*1e9, Watt2dBm(NL))
plot(lamb*1e9, Watt2dBm(NL_m), '--')
