function [Y, Df] = qpsk_freq_rec(X, FreqRec)
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

mu = FreqRec.mu;
muShift = FreqRec.muShift;

muVec = [];
muShift = [0 muShift length(X)];
for k = 1:length(mu)
    muVec = [muVec repmat(mu(k), 1, muShift(k+1)-muShift(k))];
end

thetaX = (angle((X(1, 2:end).*conj(X(1, 1:end-1))).^4)); % starts at k = 2
thetaY = (angle((X(2, 2:end).*conj(X(2, 1:end-1))).^4));

Df = zeros(1, 1+length(thetaX));
for k = 2:length(X)
    Df(k) = Df(1, k-1)*(1-muVec(k)) + muVec(k)/2*(thetaX(k-1)/(8*pi) + thetaY(k-1)/(8*pi));
end

% Df = filter([Rs*mu/(8*pi) 0], [1 -(1-mu)], theta);

% t = 0:length(Df)-1;
t = 0:length(X)-1;
Y = [X(1, :).*exp(-1j*2*pi*Df.*t);...
    X(2, :).*exp(-1j*2*pi*Df.*t)];

figure(1000), plot(56*Df)
drawnow