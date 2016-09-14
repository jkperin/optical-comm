function Pqpsk = ber_qpsk_imperfect_cpr(SNRdB, phivar)
%% Probability of error with imperfect carrier phase recovery
% Prabhu, V. K. (1976). PSK Performance with Imperfect Carrier Phase
% Recovery. IEEE Transactions on Aerospace and Electronic Systems, AES-12(2),
% 275–286. http://doi.org/10.1109/TAES.1976.308305.
% Inputs
% - SNRdB
% - phivar: phase error variance

Pqpsk = zeros(size(SNRdB));
for k = 1:length(SNRdB)
    SNR = 10^(SNRdB(k)/10);
    phibar = 0; % phase error mean
    Pqpsk(k) = 1/2*erfc(sqrt(SNR));

    tol = 1e-20;
    inc = Inf;
    l = 0;
    while abs(inc) >= tol
        Hl = sqrt(SNR)*exp(-SNR/2)/(sqrt(pi)*(2*l + 1))*(besseli(l, SNR/2)...
            + besseli(l+1, SNR/2));

        inc = (-1)^l*Hl*(1 - cos((2*l+1)*pi/4)*cos((2*l+1)*phibar)*exp(-(2*l+1)^2*phivar/2));

        Pqpsk(k) = Pqpsk(k) + inc; 
        l = l + 1;
    end
end