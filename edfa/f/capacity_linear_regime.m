function [SEnum, SEapprox] = capacity_linear_regime(E, Pump, Signal, spanAttdB, Namp, df)
%% Compute system capacity in linear regime for a particular EDF length and power loading
% Inputs:
% - E: instance of class EDF
% - Pump: instance of class Channels corresponding to pump
% - Signal: instance of class Channels corresponding to signals
% - spanAttdB: span attenuation in dB
% - Namp: number of amplifiers in the chain
% - df: frequency spacing used for computing ASE power
% Output:
% - SE: spectral efficiency i.e., capacity normalized by bandwidth

%% Numerical solution
ASEf = Channels(Signal.wavelength, 0, 'forward');
ASEb = Channels(Signal.wavelength, 0, 'backward');
[GaindB, ~, ~, Pase] = E.two_level_system(Pump, Signal, ASEf, ASEb, df, 50);

Pase = (Namp-1)*Pase;
Gain = 10.^(GaindB/10);
SNR = Gain.*Signal.P./Pase;
SEnum = (GaindB >= spanAttdB).*log2(1 + SNR); 

%% Semi-analytical solution
GaindB = E.semi_analytical_gain(Pump, Signal);
Pase = (Namp-1)*analytical_ASE_PSD(E, Pump, Signal)*df;   
Gain = 10.^(GaindB/10);
SNR = Gain.*Signal.P./Pase;
SEapprox = (GaindB >= spanAttdB).*log2(1 + SNR);   

