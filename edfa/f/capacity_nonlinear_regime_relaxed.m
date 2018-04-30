function [SE, dSE, SElamb] = capacity_nonlinear_regime_relaxed(X, E, Pump, Signal, problem)
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
% > .Gap: SNR gap to capacity in linear units
% > .step_approx: handle function to approximate step function using in
% selecting on/off channels
% > .diff_step_approx: first derivative of step approximation function
% > .excess_noise: excess noise at each wavelength 
% > .nonlinear_coeff: cell array of length 3 where D{l+2}, l = -1, 0, 1 is 
% a square matrix of the GN model nonlinear coefficients
% > .epsilon: factor to scale the nonlinear noise power from 1 span to
% "Namp" spans. 
% Output:
% - SE: negative spectral efficiency in bits/s/Hz i.e., capacity normalized by bandwidth
% - dSE: gradient of spectral efficiency with respect to X (E.L and power
% in dBm)

% Unpack parameters
spanAttdB = problem.spanAttdB;
Namp = problem.Namp;
df = problem.df;
Gap = problem.Gap;
nsp = problem.excess_noise; 
step_approx = problem.step_approx;
D = problem.nonlinear_coeff;
epsilon = problem.epsilon;

% Unpack optimization variables
E.L = X(1);
Signal.P = dBm2Watt(X(2:end));

S_and_ASE = Signal;
A = 10^(mean(spanAttdB)/10);
a = (A-1)/A;
NF = 2*a*nsp;
Acc_ASE = (Namp-1)*df*NF.*Signal.Ephoton;
S_and_ASE.P = S_and_ASE.P + Acc_ASE;

% Compute Gain using semi-analytical model
[GaindB, ~, ~, dGaindB, dGaindB_L] = E.semi_analytical_gain(Pump, S_and_ASE);

%% Nonlinear noise at the nth channel
offChs = (GaindB < spanAttdB);
P = Signal.P.*(10.^(spanAttdB/10)); % optical power after gain flattening
P(offChs) = 0; % by assigning zero power, NL noise for off channels is not computed

if nargout > 1 % gradient was requested
    [dNL, NL] = GN_model_noise_gradient(P, D);
    dNL = dNL*Namp^(1+epsilon); % Scale gradient to "Namp" spans
    
    % Multiply by gain, since dNL is originally computed with respect to 
    % the launch power, while the optimzation is done with the input power
    % to the amplifier
    dNL = dNL.*(10.^(spanAttdB/10)); 
else
    NL = GN_model_noise(P, D);
end

% Scale NL noise to "Namp" spans
NL = NL*Namp^(1+epsilon);

%% Relaxations: (i) NF is gain independent, (ii) step function approximation
SNR = Signal.P./(Namp*df*NF.*Signal.Ephoton + NL);
SElamb = 2*log2(1 + Gap*SNR).*step_approx(GaindB - spanAttdB);
SE = -sum(SElamb);

%% Computes gradient
if nargout > 1 % gradient was requested      
    % SNR gradient (validated numerically with accuray 1e-6)
    S = step_approx(GaindB - spanAttdB);
    Ninv = 1./(Namp*df*NF.*Signal.Ephoton + NL);
    dSNR = diag(Ninv) - dNL.*(Signal.P.*Ninv.^2); % SNR gradient
    % N x N matrix, where dSNR(i,j) = partial SNR_i / partial P_j
    
    diff_step_approx = problem.diff_step_approx;
        
    % SE gradient with respect to EDF length
    dSE_L = -2*sum(dGaindB_L.'.*diff_step_approx(GaindB - spanAttdB).*log2(1 + Gap*SNR));
    
    % SE gradient with respect to channel power
    %     dSElin = -2/log(2)*sum(dSNR.*(S./(1 + SNR)), 2); % Considering dG/dP = 0
    dSE_P = -2/log(2)*sum(dGaindB.*(log(1 + Gap*SNR).*diff_step_approx(GaindB - spanAttdB))  + Gap*dSNR.*(S./(1 + Gap*SNR)), 2); % Considering non-zero gain gradient
    dSE_PdBm = (log(10)/10)*(Signal.P.'.*dSE_P); % converts to derivative with power in dBm
    
    dSE = [dSE_L; dSE_PdBm];
%     dSE = dSE_PdBm;
end
