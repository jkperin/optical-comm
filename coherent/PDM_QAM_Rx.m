function [Y AGCgain] = PDM_QAM_Rx(Erec,X,Rx,t,omega)

E1rec = Erec(1,:);      % received signal in x pol.
E2rec = Erec(2,:);      % received signal in y pol.

rPBS = 10^(Rx.PolSplit.Rext/10);              % extinction ratio of PBS

% split signal + ASE
if strcmp(Rx.PolSplit.sig,'PBS')
    E11 = sqrt(rPBS/(1+rPBS))*E1rec;   % sig + ASE field in x pol. at output 1 of PBS
    E21 = sqrt(1/(1+rPBS))*E2rec;      % sig + ASE field in y pol. at output 1 of PBS
    E12 = sqrt(1/(1+rPBS))*E1rec;      % sig + ASE field in x pol. at output 2 of PBS
    E22 = sqrt(rPBS/(1+rPBS))*E2rec;   % sig + ASE field in y pol. at output 2 of PBS
elseif strcmp(Rx.PolSplit.sig,'3dB')
    E11 = sqrt(1-f3dB)*E1rec;          % sig + ASE field in x pol. at output 1 of 3-dB coupler
    E21 = sqrt(1-f3dB)*E2rec;          % sig + ASE field in y pol. at output 1 of 3-dB coupler
    E12 = sqrt(f3dB)*E1rec;            % sig + ASE field in x pol. at output 2 of 3-dB coupler
    E22 = sqrt(f3dB)*E2rec;            % sig + ASE field in y pol. at output 2 of 3-dB coupler
end

% create LO field (no phase noise)
PLO = 1e-3*10^(Rx.PLOdBm/10);                                            % Total LO power in two polarizations (W)
EL1 = sqrt(PLO/2)*ones(size(t));    % LO field in x polarization at PBS or BS input
EL2 = sqrt(PLO/2)*ones(size(t));    % LO field in y polarization at PBS or BS input

% perform polarization splitting of LO (ignore phase shifts in beamsplitters, which don't affect performance)
if strcmp(Rx.PolSplit.LO,'PBS')
    EL11 = sqrt(rPBS/(1+rPBS))*EL1;   % LO field in x pol. at output 1 of PBS
    EL21 = sqrt(1/(1+rPBS))*EL2;      % LO field in y pol. at output 1 of PBS
    EL12 = sqrt(1/(1+rPBS))*EL1;      % LO field in x pol. at output 2 of PBS
    EL22 = sqrt(rPBS/(1+rPBS))*EL2;   % LO field in y pol. at output 2 of PBS
elseif strcmp(Rx.PolSplit.LO,'3dB')
    EL11 = sqrt(1-f3dB)*EL1;          % LO field in x pol. at output 1 of 3-dB coupler
    EL21 = sqrt(1-f3dB)*EL2;          % LO field in y pol. at output 1 of 3-dB coupler
    EL12 = sqrt(f3dB)*EL1;            % LO field in x pol. at output 2 of 3-dB coupler
    EL22 = sqrt(f3dB)*EL2;            % LO field in y pol. at output 2 of 3-dB coupler
end

fS = Rx.Hybrid.fS;
fL = Rx.Hybrid.fL;
fI = Rx.Hybrid.fI;
fQ = Rx.Hybrid.fQ;

% combine signal + ASE and LO in 90-degree hybrid
phiI1 = pi/180*Rx.Hybrid.phiI01deg - omega*1e-12*Rx.Hybrid.tauIps;     % phase shift, including delay, in I branch of pol. branch 1 (rad)
phiQ1 = pi/180*Rx.Hybrid.phiQ01deg - omega*1e-12*Rx.Hybrid.tauQps;     % phase shift, including delay, in Q branch of pol. branch 1 (rad)
phiI2 = pi/180*Rx.Hybrid.phiI02deg - omega*1e-12*Rx.Hybrid.tauIps;     % phase shift, including delay, in I branch of pol. branch 2 (rad)
phiQ2 = pi/180*Rx.Hybrid.phiQ02deg - omega*1e-12*Rx.Hybrid.tauQps;     % phase shift, including delay, in Q branch of pol. branch 2 (rad)
% polarization branch 1, I
E11Ia = sqrt(1-fS)*sqrt(1-fI)*E11 + 1i*sqrt(1-fL)*sqrt(fI)*ifft(exp(1i*phiI1).*fft(EL11)); % field in x polarization at first I output in pol. branch 1
E21Ia = sqrt(1-fS)*sqrt(1-fI)*E21 + 1i*sqrt(1-fL)*sqrt(fI)*ifft(exp(1i*phiI1).*fft(EL21)); % field in y polarization at first I output in pol. branch 1
E11Ib = 1i*sqrt(1-fS)*sqrt(fI)*E11 + sqrt(1-fL)*sqrt(1-fI)*ifft(exp(1i*phiI1).*fft(EL11)); % field in x polarization at second I output in pol. branch 1
E21Ib = 1i*sqrt(1-fS)*sqrt(fI)*E21 + sqrt(1-fL)*sqrt(1-fI)*ifft(exp(1i*phiI1).*fft(EL21)); % field in y polarization at second I output in pol. branch 1
% polarization branch 1, Q
E11Qa = 1i*sqrt(fS)*sqrt(1-fQ)*ifft(exp(1i*phiQ1).*fft(E11)) - sqrt(fL)*sqrt(fQ)*EL11;     % field in x polarization at first Q output in pol. branch 1
E21Qa = 1i*sqrt(fS)*sqrt(1-fQ)*ifft(exp(1i*phiQ1).*fft(E21)) - sqrt(fL)*sqrt(fQ)*EL21;     % field in y polarization at first Q output in pol. branch 1
E11Qb = -sqrt(fS)*sqrt(fQ)*ifft(exp(1i*phiQ1).*fft(E11)) + 1i*sqrt(fL)*sqrt(1-fQ)*EL11;    % field in x polarization at second Q output in pol. branch 1
E21Qb = -sqrt(fS)*sqrt(fQ)*ifft(exp(1i*phiQ1).*fft(E21)) + 1i*sqrt(fL)*sqrt(1-fQ)*EL21;    % field in y polarization at second Q output in pol. branch 1
% polarization branch 2, I
E12Ia = sqrt(1-fS)*sqrt(1-fI)*E12 + 1i*sqrt(1-fL)*sqrt(fI)*ifft(exp(1i*phiI2).*fft(EL12)); % field in x polarization at first I output in pol. branch 2
E22Ia = sqrt(1-fS)*sqrt(1-fI)*E22 + 1i*sqrt(1-fL)*sqrt(fI)*ifft(exp(1i*phiI2).*fft(EL22)); % field in y polarization at first I output in pol. branch 2
E12Ib = 1i*sqrt(1-fS)*sqrt(fI)*E12 + sqrt(1-fL)*sqrt(1-fI)*ifft(exp(1i*phiI2).*fft(EL12)); % field in x polarization at second I output in pol. branch 2
E22Ib = 1i*sqrt(1-fS)*sqrt(fI)*E22 + sqrt(1-fL)*sqrt(1-fI)*ifft(exp(1i*phiI2).*fft(EL22)); % field in y polarization at second I output in pol. branch 2
% polarization branch 2, Q
E12Qa = 1i*sqrt(fS)*sqrt(1-fQ)*ifft(exp(1i*phiQ2).*fft(E12)) - sqrt(fL)*sqrt(fQ)*EL12;     % field in x polarization at first Q output in pol. branch 2
E22Qa = 1i*sqrt(fS)*sqrt(1-fQ)*ifft(exp(1i*phiQ2).*fft(E22)) - sqrt(fL)*sqrt(fQ)*EL22;     % field in y polarization at first Q output in pol. branch 2
E12Qb = -sqrt(fS)*sqrt(fQ)*ifft(exp(1i*phiQ2).*fft(E12)) + 1i*sqrt(fL)*sqrt(1-fQ)*EL12;    % field in x polarization at second Q output in pol. branch 2
E22Qb = -sqrt(fS)*sqrt(fQ)*ifft(exp(1i*phiQ2).*fft(E22)) + 1i*sqrt(fL)*sqrt(1-fQ)*EL22;    % field in y polarization at second Q output in pol. branch 2

% perform photodetection
R = Rx.PD.R;
Y1idet = R*((abs(E11Ia).^2 + abs(E21Ia).^2) - (abs(E11Ib).^2 + abs(E21Ib).^2));  % I output, pol. branch 1
Y1qdet = -R*((abs(E11Qa).^2 + abs(E21Qa).^2) - (abs(E11Qb).^2 + abs(E21Qb).^2));  % Q output, pol. branch 1
Y2idet = R*((abs(E12Ia).^2 + abs(E22Ia).^2) - (abs(E12Ib).^2 + abs(E22Ib).^2));  % I output, pol. branch 2
Y2qdet = -R*((abs(E12Qa).^2 + abs(E22Qa).^2) - (abs(E12Qb).^2 + abs(E22Qb).^2));  % Q output, pol. branch 2

% Filter by lowpass filter
if strcmp(Rx.LPF.Type,'Bessel')
    besfact = [1.000, 1.272, 1.405, 1.515, 1.621];
    [num,denom] = besself(Rx.LPF.Ord,2*pi*besfact(Rx.LPF.Ord)*Rx.LPF.BW);
    Hrx = polyval(num, 1i*omega)./polyval(denom, 1i*omega);   % case of Bessel LPF
elseif strcmp(Rx.LPF.Type,'Butter')
    [num,denom] = butter(Rx.LPF.Ord,2*pi*Rx.LPF.BW,'s');
    Hrx = polyval(num, 1i*omega)./polyval(denom, 1i*omega);   % case of Butterworth LPF
elseif strcmp(Rx.LPF.Type,'Cheby');% Chebyshev TypeII 
    R  = 20;
    Wo = cosh(1/Rx.LPF.Ord*acosh(sqrt(10^(R/10)-1)));
    [num, denom] = cheby2(Rx.LPF.Ord,R,Wo,'s');
    Hrx = polyval(num, 1i*omega/(2*pi*Rx.LPF.BW))./polyval(denom, 1i*omega/(2*pi*Rx.LPF.BW));   % case of Chebyshev type II LPF
end

deltaomega = omega(2)-omega(1);
gdrx = -diff(unwrap(phase(Hrx)))/deltaomega;                                % group delay at zero frequency
Hrxnd = Hrx.*exp(1i*omega*gdrx(1));                                         % freq resp with delay removed
Y1ifilt = real(ifft(fft(Y1idet).*Hrxnd));                                   % Inphase in x pol. at output of lowpass filter
Y1qfilt = real(ifft(fft(Y1qdet).*Hrxnd));                                   % Quad in x pol. at output of lowpass filter
Y2ifilt = real(ifft(fft(Y2idet).*Hrxnd));                                   % Inphase in y pol. at output of lowpass filter
Y2qfilt = real(ifft(fft(Y2qdet).*Hrxnd));                                   % Quad in x pol. at output of lowpass filter

x1i = real(X(1,:));
x1q = imag(X(1,:));
x2i = real(X(2,:));
x2q = imag(X(2,:));

% Scale received signal (analogous to AGC). The scaling makes the continuous-time
% signal + noise have the same variances as the discrete symbols transmitted. 


%### FIXME #### for 16QAM ###3
agcnrz = 0.900;                                                             % empirical AGC scale factor to achieve correct signal levels for NRZ
AGCgain = agcnrz*sqrt((var(x1i)+var(x1q)+var(x2i)+var(x2q))/(var(Y1ifilt)+var(Y1qfilt)+var(Y2ifilt)+var(Y2qfilt)));
Y = AGCgain*[Y1ifilt + 1i*Y1qfilt; Y2ifilt + 1i*Y2qfilt];
end
