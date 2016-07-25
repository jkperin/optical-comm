function [Y, Df] = qpsk_freq_rec(X, FreqRec,  verbose)
%% Frequency recovery algorithm for QPSK
% Based on: 
% 1. Sebastian Hoffmann, Suhas Bhandare, Timo Pfau, Olaf Adamczyk, Christian
% Wordehoff,  Ralf Peveling, Mario Porrmann, and Reinhold Noe. Frequency
% and phase estimation for coherent QPSK transmission with unlocked
% DFB lasers. IEEE Photonics Technology Letters, 20(18):1569{1571, 2008.
% 2. Savory, S. J. (2010). Digital Coherent Optical Receivers: 
% Algorithms and Subsystems. IEEE Journal of Selected Topics in Quantum 
% Electronics, 16(5), 1164–1179.
% Inputs:
% X : input signal in two polarizations [2 x N]
% FreqRec : FreqRec.{mu : adaptation speed (can be a vector for gear
% shifting), muShift : transition points for each change mu}
% Algorithm:
% Algorithm proceeds with adaptation rate mu until Ntrain. After that point,
% adaptation stops and the last value with frequency offset is used for the
% rest of the sequence

Ntrain = FreqRec.Ntrain;
Rs = FreqRec.Rs;
mu = FreqRec.mu;

Df = zeros(1, length(X));
theta = zeros(size(X));
Y = zeros(size(X));
for k = 2:length(X)
    theta(:, k) = angle((X(:, k).*conj(X(:, k-1))).^4); 
    
    if k < Ntrain
        Df(k) = Df(k-1)*(1-mu) + mu/(16*pi)*(theta(1, k) + theta(2, k));
    else
        Df(k) = Df(k -1);
    end
    
    Y(:, k) = X(:, k)*exp(-1j*2*pi*Df(k)*k);
end

if exist('verbose', 'var') && verbose
    figure(201), hold on, box on
    plot(Rs*Df/1e9)
    xlabel('Symbols')
    ylabel('Frequency offset (GHz)')
    drawnow
end