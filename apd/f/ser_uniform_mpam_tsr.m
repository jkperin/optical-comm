% Symbol error rate for M-PAM with uniform level spacing under thermal,
% shot and RIN noise (tsr).

% Input:
% mpam : signal struct (bw = bandwidth, M = number of levels)
% tx : transmitter struct (Prec = average received power, RIN = RIN in dB/Hz)
% rx : receiver struct (Gapd = apd gain, ka = avalanche efficiency, R = 
% responsivity)
% Gapd (optimal) : for running optimization it might be necessary to
% provide a different gain from the one in the rx struct. The excess noise
% parameter (Fa) is recalculated if Gapd is passed as a parameter.
% Prec (optimal) : for running optimization it might be necessary to
% provide a different average received power (Prec) from the one in the 
% tx struct. 

% Output:
% ser = symbol error ratio

function ser = ser_uniform_mpam_tsr(mpam, tx, rx, Gapd, Prec)

% Electron charge
q = 1.60217657e-19;

% If Gapd and Prec were passed as parameters, overwrite value in structs
% This is intended to be used during gain optimization
if nargin == 5
    rx.Gapd = Gapd;
    tx.Prec = Prec;
end   

rx.Fa = rx.ka*rx.Gapd + (1 - rx.ka)*(2 - 1/rx.Gapd);

% Define parameters so that total noise variance can be written in the form
% sigma^2 = alpha + beta*P + gamma*P^2

% Thermal noise parameter
alpha = rx.N0*mpam.bw;

% Shot noise parameter
beta = 0*rx.Gapd^2*rx.Fa*2*q*rx.R*mpam.bw;

% RIN parameter
gamma = 0*rx.Gapd^2*rx.R^2*2*10^(tx.RIN/10)*mpam.bw;

% Calculate power of the first non-zero level
P1 = 2*tx.Prec/(mpam.M-1);

% Vector of total noise variance in each level
k = 0:mpam.M-1;
sigma = sqrt(alpha + beta*k*P1 + gamma*k.^2*P1^2);

% Symbol error probability for first and last level (one way of error)
ser = 1/mpam.M*sum(qfunc(rx.Gapd*rx.R*P1./(2*sigma([1 end]))));

% Add symbol error probability of inner levels
ser = ser + 2/mpam.M*sum(qfunc(rx.Gapd*rx.R*P1./(2*sigma(2:end-1))));

