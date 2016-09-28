function [varphiE, nPN, nAWGN] = phase_error_variance(csi, wn, Ncpr, Delay, totalLinewidth, SNRdB, Rs, verbose)
%% Calculates PLL phase error variance using small-signal approximation
% Calculations are only valid for QPSK
% Calculations assume that receiver is LO shot noise limited
% Inputs:
% - csi: loop filter damping factor
% - wn: loop filter bandwidth
% Ncpr: number of polarizations used in CPR
% Delay: total loop delay
% totalLinewidth: total linewidth i.e., transmitter laser + LO
% SNRdB: SNR in dB
% Rs: symbol rate

SNR = 10^(SNRdB/10);
Mct = 1; % oversampling ratio to emulate continuous time
fs = Mct*Rs;
w = 2*pi*linspace(-fs/2, fs/2, 2^14);

nPN = zeros(size(wn));
nAWGN = zeros(size(wn));
varphiE = zeros(size(wn));
for k = 1:length(wn)
    numFs = [2*csi*wn(k) wn(k)^2];
    denFs = [1 0]; % removes integrator due to oscillator

    Fw = polyval(numFs, 1j*w)./polyval(denFs, 1j*w);
    H1 = 1j*w + exp(-1j*w*Delay).*Fw;    

    nPN(k) = totalLinewidth*trapz(w, 1./abs(H1).^2); % phase noise contribution
    nAWGN(k) = 1/(2*pi*Ncpr*2*SNR*Rs)*trapz(w, abs(Fw./H1).^2); % AWGN contribution
    varphiE(k) = nPN(k) + nAWGN(k); % phase error variance

    if exist('verbose', 'var') && verbose
        fprintf('Contribution of phase noise vs AWGN on PLL phase error at SNRdB = %.2f:\nPN/AWGN = %.3f\n', SNRdB,...
            nPN(k)/nAWGN(k));
    end
end
