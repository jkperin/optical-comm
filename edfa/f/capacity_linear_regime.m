function [num, approx] = capacity_linear_regime(E, Pump, Signal, spanAttdB, Namp, df)
%% Compute system capacity in linear regime for a particular EDF length and power loading
% Inputs:
% - E: instance of class EDF
% - Pump: instance of class Channels corresponding to pump
% - Signal: instance of class Channels corresponding to signals
% - spanAttdB: span attenuation in dB
% - Namp: number of amplifiers in the chain
% - df: frequency spacing used for computing ASE power
% Output:
% - num: struct containing the spectral efficiency, gain in dB, ASE power
% in W, and SNR in dB calculated using the numerical method
% - approx: struct similar to num, but calculated using semi-analytical
% methods

% Set channels with zero power with small power for gain/noise calculation
offChs = (Signal.P == 0);
Signal.P(offChs) = eps; 

%% Numerical solution
ASEf = Channels(Signal.wavelength, 0, 'forward');
ASEb = Channels(Signal.wavelength, 0, 'backward');
[GaindB, ~, ~, Pase] = E.two_level_system(Pump, Signal, ASEf, ASEb, df, 50);

Pase = (Namp-1)*Pase;
Gain = 10.^(GaindB/10);
SNR = Gain.*Signal.P./Pase;
SEnum = 2*(GaindB >= spanAttdB).*log2(1 + SNR); 

num.SE = SEnum;
num.GaindB = GaindB;
num.Pase = Pase;
num.SNRdB = 10*log10(SNR);

%% Semi-analytical solution
GaindB = E.semi_analytical_gain(Pump, Signal);
Pase = (Namp-1)*analytical_ASE_PSD(E, Pump, Signal)*df;   
Gain = 10.^(GaindB/10);
SNR = Gain.*Signal.P./Pase;
SEapprox = 2*(GaindB >= spanAttdB).*log2(1 + SNR);   

approx.SE = SEapprox;
approx.GaindB = GaindB;
approx.Pase = Pase;
approx.SNRdB = 10*log10(SNR);
