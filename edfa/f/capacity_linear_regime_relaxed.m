function [SE, SElamb] = capacity_linear_regime_relaxed(X, E, Pump, Signal, problem)
%% Compute system capacity in linear regime for a particular EDF length and power loading specificied in vector X
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
% > .excess_noise: excess noise at each wavelength 
% > .step_approx: handle function to approximate step function used in
% selecting on/off channels
% Output:
% - SE: spectral efficiency in bits/s/Hz i.e., capacity normalized by bandwidth
% - SElamb: secptreal efficiency in bits/s/Hz in each channel

% Unpack parameters
spanAttdB = problem.spanAttdB;
Namp = problem.Namp;
df = problem.df;
nsp = problem.excess_noise; 
step_approx = problem.step_approx;

% Unpact optimization variables
E.L = X(1);
Signal.P = dBm2Watt(X(2:end));

% Compute Gain using semi-analytical model
GaindB = E.semi_analytical_gain(Pump, Signal);

%% Relaxations: (i) NF is gain independent, (ii) step function approximation
A = 10^(mean(spanAttdB)/10);
a = (A-1)/A;
NF = 2*a*nsp;
SNR = Signal.P./(Namp*df*NF.*Signal.Ephoton);
SElamb = 2*log2(1 + SNR).*step_approx(GaindB - spanAttdB);
SE = sum(SElamb);
