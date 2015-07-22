clear, clc, close all

addpath ../mpam/
addpath ../f % general functions
addpath f

sim.Mct = 15;
% M-PAM
% M, Rb, leve_spacing, pshape
mpam = PAM(4, 100e9, 'equally-spaced', @(n) double(n >= 0 & n < sim.Mct));

sim.fs = mpam.Rs*sim.Mct;

%% Receiver
rx.N0 = (20e-12).^2; % thermal noise psd
rx.Id = 10e-9; % dark current
rx.R = 1; % responsivity
% Electric Lowpass Filter
% rx.elefilt = design_filter('bessel', 5, 0.5*mpam.Rs/(sim.fs/2));
rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);
% Optical Bandpass Filter
rx.optfilt = design_filter('fbg', 0, 200e9/(sim.fs/2));

soaG = soa(20, 9, 1310e-9, 20);

tx.RIN = -140;
tx.rexdB = -10;

Df = rx.elefilt.noisebw(sim.fs)/2;
Dfopt = rx.optfilt.noisebw(sim.fs);

varT = rx.N0*Df;

varShot = @(P) 2*1.60217657e-19*(rx.R*P + rx.Id)*Df;

varRIN = @(P) 10^(tx.RIN/10)*P.^2*Df; 

varSigSpont = @(P, G) 4*rx.R^2*P.*(G - 1).*10^(soaG.Fn/10)/2*1.5164e-19*Df;

varSpontSpont = @(G) 2*rx.R^2*((G - 1)*10^(soaG.Fn/10)/2*1.5164e-19).^2*Df*Dfopt;

noise_std = @(P, G) sqrt(varT + varShot(P) + varRIN(P) + varSigSpont(P, G) + varSpontSpont(G));
% noise_std = @(P, G) sqrt(varT + varSigSpont(P, G) + varSpontSpont(G));

tx.Ptx = 1e-3*10^(-21/10);
rex = 10^(tx.rexdB/10);
DP = tx.Ptx/(mpam.M-1)*(1 - rex)/(1 + rex);
P0 = 2*tx.Ptx*rex/(1 + rex);
P = (0:mpam.M-1)*DP + P0;
GdB = 1:20;
G =  10.^(GdB/10);
n = [1 2 2 1];

figure, hold on
pe = 0;
for k = 1:length(P)
    plot(GdB, log10(n(k)/mpam.M*qfunc(G*DP./noise_std(P(k)*G, G))))
    pe = pe + n(k)/mpam.M*qfunc(G*DP./noise_std(P(k)*G, G));
end
plot(GdB, log10(pe), 'k')
axis([GdB([1 end]) -10 0])

% figure
% plot(GdB, log10(qfunc(DP./noise_std(tx.Ptx, 10.^(GdB/10)))))