function [num, approx] = capacity_nonlinear_regime(E, Pump, Signal, problem)
%% Compute system capacity in linear regime for a particular EDF length and power loading
% Inputs:
% - E: instance of class EDF
% - Pump: instance of class Channels corresponding to pump
% - Signal: instance of class Channels corresponding to signals
% - problem: struct containing parameters from particular problem
% > .spanAttdB: span attenuation in dB at each signal wavelength
% > .Namp: number of amplifiers in the chain
% > .df: channel spacing. Use to compute noise power
% > .Gap: SNR gap to capacity in linear units
% > .excess_noise: excess noise at each wavelength 
% > .nonlinear_coeff: cell array of length 3 where D{l+2}, l = -1, 0, 1 is 
% a square matrix of the GN model nonlinear coefficients
% > .epsilon: factor to scale the nonlinear noise power from 1 span to
% "Namp" spans. 
% Output:
% - num: struct containing the spectral efficiency, gain in dB, ASE power
% in W, and SNR in dB calculated using the numerical method
% - approx: struct similar to num, but calculated using semi-analytical
% methods

% Unpack parameters
spanAttdB = problem.spanAttdB;
Namp = problem.Namp;
df = problem.df;
Gap = problem.Gap;
nsp = problem.excess_noise; 
D = problem.nonlinear_coeff;
epsilon = problem.epsilon;

% Set channels with zero power with small power for gain/noise calculation
offChs = (Signal.P == 0);
Signal.P(offChs) = eps; 

%% Numerical solution
ASEf = Channels(Signal.wavelength, 0, 'forward');
ASEb = Channels(Signal.wavelength, 0, 'backward');

S_and_ASE = Signal;
A = 10^(mean(spanAttdB)/10);
a = (A-1)/A;
NF = 2*a*nsp;
Acc_ASE = (Namp-1)*df*NF.*Signal.Ephoton;
S_and_ASE.P = S_and_ASE.P + Acc_ASE;
[GaindB, ~, ~, Pase] = E.propagate(Pump, S_and_ASE, ASEf, ASEb, df, 'two-level');

Gain = 10.^(GaindB/10);
Pase = Namp*Pase./Gain;

%% Nonlinear noise at the nth channel
offChs = (GaindB < spanAttdB);
P = Signal.P.*(10.^(spanAttdB/10)); % optical power after gain flattening
P(offChs) = 0;
NL = GN_model_noise(P, D);

% Scale NL noise to "Namp" spans
NL = (10.^(-spanAttdB/10)).*NL*(Namp^(1+epsilon));

% Compute capacity
SNR = Signal.P./(Pase + NL);
SEnum = 2*(GaindB >= spanAttdB).*log2(1 + Gap*SNR); 

%
num.SE = SEnum;
num.GaindB = GaindB;
num.Pase = Pase;
num.SNRdB = 10*log10(SNR);
num.NL = NL;

%% Semi-analytical solution
GaindB = E.semi_analytical_gain(Pump, S_and_ASE);
Gain = 10.^(GaindB/10);
Pase = Namp*df*(2*nsp.*(Gain-1).*Signal.Ephoton)./Gain;
Pase = max(Pase, 0); % non-negativity constraint           

SNR = Signal.P./(Pase + NL);
SEapprox = 2*(GaindB >= spanAttdB).*log2(1 + Gap*SNR);   

approx.SE = SEapprox;
approx.GaindB = GaindB;
approx.Pase = Pase;
approx.SNRdB = 10*log10(SNR);
approx.NL = NL;

