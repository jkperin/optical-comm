% BER for M-PAM with uniform level spacing under thermal, shot and RIN noie

% Input:
% mpam : signal struct (bw = bandwidth, M = number of levels)
% tx : transmitter struct (Prec = average received power, RIN = RIN in dB/Hz)
% rx : receiver struct (Gapd = apd gain, ka = avalanche efficiency, R = 
% responsivity)
% sim : simulation struct (BERtarget = target BER, rin = on/off, shot
% = on/off, levelSpacing = uniform or nonuniform)
% Gapd (optimal) : for running optimization it might be necessary to
% provide a different gain from the one in the rx struct.
% Prec (optimal) : for running optimization it might be necessary to
% provide a different average received power (Prec) from the one in the 
% tx struct. 

% Output:
% ber = bit error ratio (assuming that one symbol error leads to
% approximately one bit error)

function ber = ber_mpam(mpam, tx, rx, sim, Gapd, Prec)

% Electron charge
q = 1.60217657e-19;

% If Gapd and Prec were passed as parameters, overwrite value in structs
% This is intended to be used during gain optimization
if nargin >= 5
    rx.Gapd = Gapd;
    tx.Prec = Prec;
end   

% Define parameters so that total noise variance can be written in the form
% sigma^2 = alpha + beta*P + gamma*P^2

% Thermal noise parameter
alpha = rx.N0*mpam.bw;

% Shot noise parameter
beta = 0;
if strcmp(sim.shot, 'on')
    beta = rx.Gapd^2*rx.Fa(rx.ka, rx.Gapd)*2*q*rx.R*mpam.bw;
end

% RIN parameter
gamma = 0;
if strcmp(sim.rin, 'on')
    gamma = rx.Gapd^2*rx.R^2*2*10^(tx.RIN/10)*mpam.bw;
end

if strcmp(sim.levelSpacing, 'uniform')
    % Calculate power of the first non-zero level
    P1 = 2*tx.Prec/(mpam.M-1);   

    % Vector of total noise variance in each level
    Pk = (0:mpam.M-1)*P1;  % Power of the k-th level
    sigma = sqrt(alpha + beta*Pk + gamma*Pk.^2);

    % Symbol error probability for first and last level (one way of error)
    ser = 1/mpam.M*sum(qfunc(rx.Gapd*rx.R*P1./(2*sigma([1 end]))));

    % Add symbol error probability of inner levels
    ser = ser + 2/mpam.M*sum(qfunc(rx.Gapd*rx.R*P1./(2*sigma(2:end-1))));

elseif strcmp(sim.levelSpacing, 'nonuniform')
    % Scale power to the desired transmitted power. This remains the
    % optimal level spacing only when tx.Prec == mean(a)
    [a, b] = calcOptLevelSpacing(mpam, tx, rx, sim);
    
    Pk = a/mean(a)*tx.Prec;  % Power of the k-th level
    bk = b/mean(a)*tx.Prec;  % new decision thresholds

    % append an element to make indexing easier
    bk = [NaN; bk];

    % Variance at different power levels
    sigma = sqrt(alpha + beta*Pk + gamma*Pk.^2);

    % Symbol error probability for first and last level (one way of error)
    ser = 1/mpam.M*sum(qfunc(rx.Gapd*rx.R*abs(Pk([1 end]) - bk([2 end]))./(sigma([1 end]))));

    % Add symbol error probability of inner levels
    ser = ser + 1/mpam.M*sum(qfunc(rx.Gapd*rx.R*(Pk(2:end-1) - bk(2:end-1))./(sigma(2:end-1))) + ...
                       qfunc(rx.Gapd*rx.R*(bk(3:end) - Pk(2:end-1))./(sigma(2:end-1))));
else 
    error('Invalid Option!');
end

% Assuming Gray coding, one symbol error leads to approximately one bit
% error. This approximation might not be accurate for uniform level spacing
% since the noise at higher levels is significantly higher. THis
% approximation also breaks at small SNRs.
ber = ser/log2(mpam.M);

