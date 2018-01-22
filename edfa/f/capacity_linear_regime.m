function [num, approx] = capacity_linear_regime(E, Pump, Signal, problem)
%% Compute system capacity in linear regime for a particular EDF length and power loading
% Inputs:
% - E: instance of class EDF
% - Pump: instance of class Channels corresponding to pump
% - Signal: instance of class Channels corresponding to signals
% - problem: struct containing parameters from particular problem
% > .spanAttdB: span attenuation in dB at each signal wavelength
% > .Namp: number of amplifiers in the chain
% > .df: channel spacing. Use to compute noise power
% > .excess_noise: excess noise at each wavelength. This is only used for
% the semi-analytical method
% Output:
% - num: struct containing the spectral efficiency, gain in dB, ASE power
% in W, and SNR in dB calculated using the numerical method
% - approx: struct similar to num, but calculated using semi-analytical
% methods

% Unpack parameters
spanAttdB = problem.spanAttdB;
Namp = problem.Namp;
df = problem.df;
nsp = problem.excess_noise; 

% Set channels with zero power with small power for gain/noise calculation
offChs = (Signal.P == 0);
Signal.P(offChs) = eps; 

%% Numerical solution
ASEf = Channels(Signal.wavelength, 0, 'forward');
ASEb = Channels(Signal.wavelength, 0, 'backward');
[GaindB, ~, ~, Pase] = E.two_level_system(Pump, Signal, ASEf, ASEb, df);

Gain = 10.^(GaindB/10);
Pase = Namp*Pase./Gain;
SNR = Signal.P./Pase;
SEnum = 2*(GaindB >= spanAttdB).*log2(1 + SNR); 

num.SE = SEnum;
num.GaindB = GaindB;
num.Pase = Pase;
num.SNRdB = 10*log10(SNR);

%% Semi-analytical solution
GaindB = E.semi_analytical_gain(Pump, Signal);
Gain = 10.^(GaindB/10);
Pase = Namp*df*(2*nsp.*((Gain-1)./Gain).*Signal.Ephoton);
Pase = max(Pase, 0); % non-negativity constraint           

SNR = Signal.P./Pase;
SEapprox = 2*(GaindB >= spanAttdB).*log2(1 + SNR);   

approx.SE = SEapprox;
approx.GaindB = GaindB;
approx.Pase = Pase;
approx.SNRdB = 10*log10(SNR);
