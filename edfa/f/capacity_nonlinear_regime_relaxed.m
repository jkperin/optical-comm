function [SE, SElamb] = capacity_nonlinear_regime_relaxed(X, E, Pump, Signal, problem)
%% Compute system capacity in nonlinear regime for a particular EDF length and power loading specificied in vector X
% X(1) is the EDF length, and X(2:end) has the signal power in dBm at each
% wavelength. Simulations assume ideal gain flatenning, resulting in the
% simplified capacity formula
% Inputs:
% - X: vector containing EDF length and power load
% - E: instance of class EDF
% - Pump: instance of class Channels corresponding to pump
% - Signal: instance of class Channels corresponding to signals
% - problem: struct containing parameters from particular problem
% > .spanAttdB: span attenuation in dB at each signal wavelength
% > .Namp: number of amplifiers in the chain
% > .df: channel spacing. Use to compute noise power
% > .step_approx: handle function to approximate step function using in
% selecting on/off channels
% > .excess_noise: excess noise at each wavelength 
% > .nonlinear_coeff: cell array of length 3 where D{l+2}, l = -1, 0, 1 is 
% a square matrix of the GN model nonlinear coefficients
% > .epsilon: factor to scale the nonlinear noise power from 1 span to
% "Namp" spans. 
% Output:
% - SE: spectral efficiency in bits/s/Hz i.e., capacity normalized by bandwidth
% - SElamb: spectral efficiency in bits/s/Hz in each wavelength 

% Unpack parameters
spanAttdB = problem.spanAttdB;
Namp = problem.Namp;
df = problem.df;
nsp = problem.excess_noise; 
step_approx = problem.step_approx;
D = problem.nonlinear_coeff;
epsilon = problem.epsilon;

% Unpack optimization variables
E.L = X(1);
Signal.P = dBm2Watt(X(2:end));

% Compute Gain using semi-analytical model
GaindB = E.semi_analytical_gain(Pump, Signal);

%% Nonlinear noise at the nth channel
offChs = (GaindB < spanAttdB);
P = Signal.P.*(10.^(spanAttdB/10)); % optical power after gain flattening
P(offChs) = 0;
NL = GN_model_noise(P, D);

% Scale NL noise to "Namp" spans
NL = NL*Namp^(1+epsilon);

%% Relaxations: (i) NF is gain independent, (ii) step function approximation
A = 10^(mean(spanAttdB)/10);
a = (A-1)/A;
NF = 2*a*nsp;
SNR = Signal.P./(Namp*df*NF.*Signal.Ephoton + NL);
SElamb = 2*log2(1 + SNR).*step_approx(GaindB - spanAttdB);
SE = sum(SElamb);