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
% > .nsp: excess noise. Either a fixed value or calculated analytically 
% from fiber parameters (power independent)
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
D = problem.nonlinear_coeff;

% Set channels with zero power with small power for gain/noise calculation
offChs = (Signal.P == 0);
Signal.P(offChs) = eps; 

%% Numerical solution
ASEf = Channels(Signal.wavelength, 0, 'forward');
ASEb = Channels(Signal.wavelength, 0, 'backward');
[GaindB, ~, ~, Pase] = E.two_level_system(Pump, Signal, ASEf, ASEb, df);

Gain = 10.^(GaindB/10);
Pase = Namp*Pase./Gain;

%% Nonlinear noise at the nth channel
tic
offChs = (GaindB < spanAttdB);
P = Signal.P.*(10.^(spanAttdB/10)); % optical power after gain flattening
P(offChs) = 0;
NL = zeros(1, Signal.N);
for n = 1:Signal.N
    if P(n) == 0
        continue
    end
    for n1 = 1:Signal.N
        if P(n1) == 0
            continue
        end
        for n2 = 1:Signal.N
            for l = -1:1
                idx = n1 + n2 - n + l;
                if idx < 1 || idx > Signal.N
                    continue
                end                
                if P(n2) == 0 || P(idx) == 0
                    continue
                end
                
                i = Signal.N - (n1 - n);
                j = Signal.N + (n2 - n);
                NL(n) = NL(n) + P(n1)*P(n2)*P(idx)*D{l+2}(i, j);
            end
        end
    end
end

toc

% Scale NL noise to "Namp" spans
epsilon = 0.05; % From Fig. 17 of P. Poggiolini and I. Paper, “The GN Model
% of Non-Linear Propagation in Uncompensated Coherent Optical Systems,” 
% J. Lightw. Technol., vol. 30, no. 24, pp. 3857–3879, 2012.
NL = NL*Namp^(1+epsilon);

% Compute capacity
SNR = Signal.P./(Pase + NL);
SEnum = 2*(GaindB >= spanAttdB).*log2(1 + SNR); 

%
num.SE = SEnum;
num.GaindB = GaindB;
num.Pase = Pase;
num.SNRdB = 10*log10(SNR);
num.NL = NL;

%% Semi-analytical solution
GaindB = E.semi_analytical_gain(Pump, Signal);
Gain = 10.^(GaindB/10);
Pase = Namp*df*(2*nsp.*(Gain-1).*Signal.Ephoton)./Gain;
Pase = max(Pase, 0); % non-negativity constraint           

SNR = Signal.P./(Pase + NL);
SEapprox = 2*(GaindB >= spanAttdB).*log2(1 + SNR);   

approx.SE = SEapprox;
approx.GaindB = GaindB;
approx.Pase = Pase;
approx.SNRdB = 10*log10(SNR);
approx.NL = NL;

