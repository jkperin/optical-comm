function [varphiE, nPN, nFlicker, nAWGN] = phase_error_variance_with_flicker(ka, csi, wn, Ncpr, Delay, totalLinewidth, SNRdB, Rs, verbose)
%% Calculates PLL phase error variance using small-signal approximation
%% Same as phase_error_variance.m, but including an additional term due to flicker noise
% Calculations are only valid for QPSK
% Calculations assume that receiver is LO shot noise limited
% Inputs:
% - ka: ka/f^3 is the double-sided phase noise PSD due to flicker noise
% - csi: loop filter damping factor
% - wn: loop filter bandwidth
% - Ncpr: number of polarizations used in CPR
% - Delay: total loop delay
% - totalLinewidth: total linewidth i.e., transmitter laser + LO
% - SNRdB: SNR in dB
% - Rs: symbol rate

warning('off', 'MATLAB:integral:MaxIntervalCountReached')

SNR = 10^(SNRdB/10);
Fmax = Rs; % Maximum frequency used in integration
% Note: integral from matlab supports -Inf, Inf limits, but the integration
% fails occasionally when using these limits for the integrals in this
% function

nPN = zeros(size(wn));
nAWGN = zeros(size(wn));
nFlicker = zeros(size(wn));
varphiE = zeros(size(wn));
for k = 1:length(wn)
    numFs = [2*csi*wn(k) wn(k)^2];
    denFs = [1 0]; % removes integrator due to oscillator

    Fw = @(w) polyval(numFs, 1j*w)./polyval(denFs, 1j*w);
    Hpn = @(w) 1./abs(1j*w + exp(-1j*w*Delay).*Fw(w)).^2;
    Hawgn = @(w) abs(Fw(w)./(1j*w + exp(-1j*w*Delay).*Fw(w))).^2;

    nPN(k) = totalLinewidth*integral(Hpn, -Fmax, Fmax); % phase noise contribution
    nFlicker(k) = 2*(2*pi)^2*ka*integral(@(w) Hpn(w)./w, eps, Fmax);
    nAWGN(k) = 1/(2*pi*Ncpr*2*SNR*Rs)*integral(Hawgn, -Fmax, Fmax); % AWGN contribution
    varphiE(k) = nPN(k) + nFlicker(k) + nAWGN(k); % phase error variance
       
    if exist('verbose', 'var') && verbose
        fprintf('Contribution on phase error variance at SNRdB = %.2f:\nPN/AWGN = %.3f\nFlicker/AWGN = %.3f\n', SNRdB,...
            nPN(k)/nAWGN(k), nFlicker(k)/nAWGN(k));
        
        figure(320)
        w = 2*pi*linspace(-10e9, 10e9, 2^14);
        subplot(211)
        plot(w/(2*pi*1e9), 10*log10(Hpn(w)))
        xlabel('Frequency (GHz)')
        ylabel('Amplitude response')
        title('Phase noise component integrand')
        subplot(212)
        plot(w/(2*pi*1e9), 10*log10(Hawgn(w)));
        xlabel('Frequency (GHz)')
        ylabel('Amplitude response')
        title('AWGN component integrand')
    end
end
