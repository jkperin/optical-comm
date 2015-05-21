%% Calculate Optimal level spacing for M-PAM with thermal,
%% shot and intensity noises. The levels are chosen so that the Q-factors
%% at all levels are the same (sub-optimal).

% Input:
% mpam : signal struct (bw = bandwidth, M = number of levels)
% tx : transmitter struct (RIN = RIN in dB/Hz)
% rx : receiver struct (Gapd = apd gain, ka = avalanche efficiency, Fa =
% excess noise, R = responsivity)
% sim : simulation paramters (struct) (SERtarget = target SER)
% Gapd (optimal) : for running optimization it might be necessary to
% provide a different gain from the one in the rx struct. The excess noise
% parameter (Fa) is recalculated if Gapd is passed as a parameter.

% Output:
% a : levels
% b : decision threshold

function [a, b] = calcOptLevelSpacingAmpl(mpam, tx, soa, rx, sim, Gsoa)

% Electron charge
q = 1.60217657e-19;

% If Gapd was passed as a parameter, overwrite value in struct
% This is intended to be used during gain optimization
if nargin == 6
    soa.G = Gsoa;
end   

% Define parameters so that total noise variance can be written in the form
% sigma^2 = alpha + beta*P + gamma*P^2
% Thermal noise parameter
alpha = rx.N0*mpam.bw + 4*(soa.Seq(soa.G)*mpam.Rs/2)^2;

% Shot noise parameter
beta = 2*soa.G*mpam.Rs*soa.Seq(soa.G);

% RIN parameter
gamma = 0;
if strcmp(sim.rin, 'on')
    gamma = rx.Gapd^2*rx.R^2*2*10^(tx.RIN/10)*mpam.bw;
end

% xi is the argument of the Q function. In the non-uniform level
% spacing case all the levels have the same Q-factor.
xi = qfuncinv(log2(mpam.M)*sim.BERtarget*mpam.M/(2*(mpam.M-1)));

delta = xi/(soa.G*rx.R);

% Calculate levels ak, k = 0,..., M-1
a = zeros(mpam.M, 1);
sigma = zeros(mpam.M, 1);
for k = 2:mpam.M
    sigma(k-1) = sqrt(alpha + beta*a(k-1) + gamma*a(k-1)^2);
    c1 = a(k-1) + sigma(k-1)*delta;
    
    a(k) = -a(k-1) + (2*c1 + delta^2*beta)/(1-delta^2*gamma);
end
    
% Decision threshold
b = a(1:end-1) + sigma(1:end-1)*delta;
