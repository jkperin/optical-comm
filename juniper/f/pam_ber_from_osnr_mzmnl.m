function ber = pam_ber_from_osnr_mzmnl(M, OSNRdB, noiseBW, Vset, rexdB)
%% Caculates BER from PAM transmission in IM-DD system including penalty due to MZM nonlinearity.
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
mpam = mpam.norm_levels;
mpam.b = 2*Vset(1)*(mpam.b - mean(mpam.a)) + Vset(2);
mpam.a = 2*Vset(1)*(mpam.a - mean(mpam.a)) + Vset(2);
mpamV = mpam;
mpam.a = sin(pi/2*mpam.a).^2;
mpam.b = mpam.place_thresholds();
Pmean = mean(mpam.a);

BWref = 12.5e9; % OSNR measurement bandwidth
ber = zeros(size(OSNRdB));
for k = 1:length(OSNRdB)
    OSNR = 10^(OSNRdB(k)/10);
    
    Ssp = Pmean/(2*BWref*OSNR); % Amplifier one-sided ASE PSD 
    
    % Noise std for intensity level Plevel
    noise_std = @(Plevel) sqrt(4*Plevel*Ssp*noiseBW + 0*0.65e-3);
   
    ber(k) = mpam.berAWGN(noise_std);
end

figure(233), clf, hold on, box on
t = linspace(0, 1);
plot(t, sin(pi/2*t).^2, 'k');
plot((mpamV.a*[1 1]).', [zeros(1, M); mpam.a.'], 'k');
plot([zeros(1, M); mpamV.a.'], (mpam.a*[1 1]).', 'k')
plot([0 1], (mpam.b*[1 1]).', ':k')
xlabel('Driving signal')
ylabel('Resulting power levels')
axis([0 1 0 1])
1;