function ber = pam_ber_from_osnr(M, OSNRdB, noiseBW, rexdB)
%% Caculates BER from PAM transmission in IM-DD system.
% This calculation assumes that signal-spontaneous beat noise is dominant.
% It disregards spontaneous-spontaneous beat noise and shot noise from
% local oscillator in the case of using coherent receiver.
% Inputs:
% - M: PAM order
% - OSNRdB: OSNR in dB
% - noiseBW: noise bandwidth of receiver filtering. If noiseBW =
% symbolRate/2, then noise enhacement penalty is not included
% - rexdB: extinction ratio in dB. Defined as 10*log10(Pmin/Pmax) (dB)

if not(exist('rexdB', 'var'))
    rexdB = -Inf;
end

mpam = PAM(M, 1);
mpam = mpam.adjust_levels(1, rexdB);
Pmean = mean(mpam.a);

BWref = 12.5e9; % OSNR measurement bandwidth
ber = zeros(size(OSNRdB));
for k = 1:length(OSNRdB)
    OSNR = 10^(OSNRdB(k)/10);
    
    N0 = Pmean/(BWref*OSNR); % Amplifier one-sided ASE PSD 
    
    % Noise std for intensity level Plevel
    noise_std = @(Plevel) sqrt(2*Plevel*N0*noiseBW);
   
    ber(k) = mpam.berAWGN(noise_std);
end


