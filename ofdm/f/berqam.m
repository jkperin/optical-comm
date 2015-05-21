% Calculates bit error probabilities for square (and rectangular)
% constellations as well as cross constellations. The Matlab function berawgn
% only works for rectangular/square constellations.

% [1] J. G. Smith, "Odd-bit quadrature amplitude shift keying," IEEE Trans.
% Commun., vol. COMM-23, no. 3, pp. 385–389, Mar. 1975.
% [2] P. K. Vitthaladevuni, et al. "BER Computation for Cross QAM
% Constellations," IEEE TRANSACTIONS ON WIRELESS COMMUNICATIONS, Nov. 2005.
% [3] Cho, K., and Yoon, D., "On the general BER expression of one- and 
% two-dimensional amplitude modulations", IEEE Trans. Commun., 2002.

% CS must be an integer power of 2
% SNR can be a vector

% if SNR is a vector ber is a vector of same length

function ber = berqam(CS, SNR)

if CS == 2
    % In case the power allocation sets one constellation to 2
    ber = berawgn(SNR, 'pam', CS);
        
elseif mod(log2(CS), 2) == 0 || CS == 8
    % if number of bits/symbol is even or CS == 8 (i.e., square constellation 
    % or rectangular contellation, respectively) uses matlab built-in function
    
    ber = berawgn(SNR-10*log10(log2(CS)), 'qam', CS); % based on [3]
    
else
    % otherwise we have a cross constellation. 
    % Calculates BER using the Smith approximation [1] (also see [2]) 
    
    % Gray penalty (i.e., penalty due to the imperfect Gray mapping on
    % cross constellations). 
    % Gp = 7/6 for M = 32 and
    % Gp = (1 + 1/sqrt(2*M) + 1/(3*M)), for M > 32
    Gp = @(M) 7/6*(M == 32) + (M ~= 32)*(1 + 1/sqrt(2*M) + 1/(3*M)); 
    
    % N is the average number of nearest neighbors for a symbol in the constellation
    N = @(M) 4 - 6/sqrt(2*M);

    % Smith approximation
    ber = Gp(CS)*N(CS)/log2(CS)*qfunc(sqrt(3*10.^(SNR/10)/(31/32*CS-1)));
end
    

    
