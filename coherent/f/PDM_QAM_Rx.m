function [Y, SNR_est, ACG] = PDM_QAM_Rx(Erec, M, Rx, sim)

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

% Create LO field 
ELO = Rx.LO.cw(sim);
EL1 = sqrt(1/2)*ELO;    % LO field in x polarization at PBS or BS input
EL2 = sqrt(1/2)*ELO;    % LO field in y polarization at PBS or BS input

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

% Hybrid parameters
fS = Rx.Hybrid.fS;
fL = Rx.Hybrid.fL;
fI = Rx.Hybrid.fI;
fQ = Rx.Hybrid.fQ;

% combine signal + ASE and LO in 90-degree hybrid
omega = 2*pi*sim.f;
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
Y1idet = -Rx.PD.detect([E11Ia; E21Ia], sim.fs, 'gaussian')...
    + Rx.PD.detect([E11Ib; E21Ib], sim.fs, 'gaussian'); % I output, pol. branch 1

Y1qdet = Rx.PD.detect([E11Qa; E21Qa], sim.fs, 'gaussian')...
    - Rx.PD.detect([E11Qb; E21Qb], sim.fs, 'gaussian'); % Q output, pol. branch 1

Y2idet = -Rx.PD.detect([E12Ia; E22Ia], sim.fs, 'gaussian')...
    + Rx.PD.detect([E12Ib; E22Ib], sim.fs, 'gaussian'); % I output, pol. branch 2

Y2qdet = Rx.PD.detect([E12Qa; E22Qa], sim.fs, 'gaussian')...
    - Rx.PD.detect([E12Qb; E22Qb], sim.fs, 'gaussian'); % Q output, pol. branch 2

% Add thermal noise
Y1idet = Y1idet + sqrt(Rx.N0*sim.fs/2)*randn(size(Y1idet));
Y1qdet = Y1qdet + sqrt(Rx.N0*sim.fs/2)*randn(size(Y1qdet));
Y2idet = Y2idet + sqrt(Rx.N0*sim.fs/2)*randn(size(Y2idet));
Y2qdet = Y2qdet + sqrt(Rx.N0*sim.fs/2)*randn(size(Y2qdet));

% Form output
Y = [Y1idet + 1j*Y1qdet; Y2idet + 1j*Y2qdet];

Pin = mean(abs(E11Ia).^2 + abs(E21Ia).^2);
Df = Rx.ADC.filt.noisebw(sim.fs)/2; % receiver equivalent noise bandwidth (one-sided)
varShot = 4*Rx.PD.varShot(Pin, Df); % x 4 because of 4 photodiodes are used per polarization

% Thermal noise is assumed to be per real dimension
SNR_est = 10*log10(0.5*mean(abs(Y(1, :)).^2)/(varShot + 2*Rx.N0*Df) - 1);

% AGC
% Normalize signals so that constellations be in original form
Enorm = mean(abs(qammod(0:M-1, M)).^2);
ACG = sqrt(Enorm./[mean(abs(Y(1, :)).^2), mean(abs(Y(2, :)).^2)]);

Y(1, :) = ACG(1)*Y(1, :);
Y(2, :) = ACG(1)*Y(2, :);