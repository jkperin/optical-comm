function yk = coherent_imdd_rx(Erx, mpam, TxLaser, Rx, sim)
%% Use Dual Pol coherent receiver to obtain direct detection
% PD.detect(Ein: Input electric field, fs: sampling rate of samples in Ein,
% noise statistics {'no noise', 'gaussian'}, N0: one-sided PSD of thermal
% noise)
LO = laser(TxLaser.wavelength, Rx.PlodBm, -150, 200e3);

PolSplit.sig = 'PBS';
PolSplit.LO = '3dB';

rPBS = 10^(-30/10);              % extinction ratio of PBS

% split signal + ASE
if strcmp(PolSplit.sig,'PBS')
    E1rec = Erx(1,:);      % received signal in x pol.
    E2rec = Erx(2,:);      % received signal in y pol.
    Es1 = [sqrt(rPBS/(1+rPBS))*E1rec;   % sig + ASE field in x pol. at output 1 of PBS
           sqrt(1/(1+rPBS))*E2rec];      % sig + ASE field in y pol. at output 1 of PBS
    Es2 = [sqrt(1/(1+rPBS))*E1rec;      % sig + ASE field in x pol. at output 2 of PBS   
           sqrt(rPBS/(1+rPBS))*E2rec];   % sig + ASE field in y pol. at output 2 of PBS
elseif strcmp(PolSplit.sig,'3dB')
    [Es1, Es2] = powerSplitter(Erx, 0, 0.5, 0);
end

% Create LO field 
ELO = LO.cw(sim); % generates continuous-wave electric field in 1 pol with intensity and phase noise
ELO = [sqrt(1/2)*ELO;    % LO field in x polarization at PBS or BS input
       sqrt(1/2)*ELO];    % LO field in y polarization at PBS or BS input

% perform polarization splitting of LO (ignore phase shifts in beamsplitters, which don't affect performance)
if strcmp(PolSplit.LO,'PBS')
    EL1 = [sqrt(rPBS/(1+rPBS))*ELO(1, :);   % LO field in x pol. at output 1 of PBS
           sqrt(1/(1+rPBS))*ELO(2, :)];      % LO field in y pol. at output 1 of PBS
    EL2 = [sqrt(1/(1+rPBS))*ELO(1, :);      % LO field in x pol. at output 2 of PBS
            sqrt(rPBS/(1+rPBS))*ELO(2, :)];   % LO field in y pol. at output 2 of PBS
elseif strcmp(PolSplit.LO,'3dB')
    [EL1, EL2] = powerSplitter(ELO, 0, 0.5, 0);
end

% 90deg Hybrid
% Note: fields are 2 x N matrices, where each row is one polarization
% Power splitter of signal and LO for pol 1
[Es1i, Es1q] = powerSplitter(Es1, 0, 0.5, 0);
[EL1i, EL1q] = powerSplitter(EL1, 0, 0.5, 0);
% Power splitter of signal and LO for pol 2
[Es2i, Es2q] = powerSplitter(Es2, 0, 0.5, 0);
[EL2i, EL2q] = powerSplitter(EL2, 0, 0.5, 0);

% Note: last indices 1 and 2 refer to output ports rather than polarizations
% Power splitter of I of pol 1
[E1i1, E1i2] = powerSplitter(Es1i, EL1i*exp(-1j*pi/2), 0.5, 0);
% Power splitter of Q of pol 1
[E1q1, E1q2] = powerSplitter(Es1q, EL1q, 0.5, 0);

% Power splitter of I of pol 2
[E2i1, E2i2] = powerSplitter(Es2i, EL2i*exp(-1j*pi/2), 0.5, 0);
% Power splitter of Q of pol 2
[E2q1, E2q2] = powerSplitter(Es2q, EL2q, 0.5, 0);

% perform photodetection
I(1, :) = Rx.PD.detect(E1i1, sim.fs, 'gaussian')...
    - Rx.PD.detect(E1i2, sim.fs, 'gaussian'); % I output, pol. branch 1

I(2, :) = Rx.PD.detect(E1q1, sim.fs, 'gaussian')...
    - Rx.PD.detect(E1q2, sim.fs, 'gaussian'); % Q output, pol. branch 1

I(3, :) = Rx.PD.detect(E2i1, sim.fs, 'gaussian')...
    - Rx.PD.detect(E2i2, sim.fs, 'gaussian'); % I output, pol. branch 2

I(4, :) = Rx.PD.detect(E2q1, sim.fs, 'gaussian')...
    - Rx.PD.detect(E2q2, sim.fs, 'gaussian'); % Q output, pol. branch 2

% Add thermal noise
I = I + sqrt(Rx.N0*sim.fs/2)*randn(size(I));

%% ADC
% ADC performs filtering, quantization, and downsampling
% For an ideal ADC, ADC.ENOB = Inf
% Rx.ADC.offset = -1;
Ixik = adc(I(1,:), Rx.ADC, sim);
Ixqk = adc(I(2,:), Rx.ADC, sim);
Iyik = adc(I(3,:), Rx.ADC, sim);
Iyqk = adc(I(4,:), Rx.ADC, sim);

%% Prefiltering. Same as done in process_pam_waveforms.
fs = Rx.ADC.fs;
f = freq_time(length(Ixik), fs);
% Filt = design_filter('butter', 5, mpam.Rs/(fs/2)); % half of the sampling rate
% 
% Hprefilt = ifftshift(Filt.H(f/fs));
% Ixik = real(ifft(fft(Ixik).*Hprefilt));
% Ixqk = real(ifft(fft(Ixqk).*Hprefilt));
% Iyik = real(ifft(fft(Iyik).*Hprefilt));
% Iyqk = real(ifft(fft(Iyqk).*Hprefilt));

%% Combine I and Q to obtain intensity
yt = abs(Ixik + 1j*Ixqk).^2 + abs(Iyik + 1j*Iyqk).^2;

% Gan control
yt = yt - mean(yt);
yt = yt*sqrt(mean(abs(mpam.a).^2)/(sqrt(2)*mean(abs(yt).^2)));

%% Antialiasing filter and reduce noise
fBW = Rx.Filt.fcnorm*sim.fs/2;
Filt = design_filter(Rx.Filt.type, Rx.Filt.order, fBW/(fs/2)); % half of the sampling rate
yt = real(ifft(fft(yt).*ifftshift(Filt.H(f/fs))));

%% Resample
% Resample in order to have integer number of samples per symbol specified by ros
[p, q] = rat(sim.ros.rxDSP*mpam.Rs/fs);
yk = resample(yt, p, q); % resample to match symbol rate at the DAC

yk = yk - mean(yk);