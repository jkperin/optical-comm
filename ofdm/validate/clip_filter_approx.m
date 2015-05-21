% test whether reversing the order of filter and clipping is a good approximation

clc, clear, close all

% initiate default random number generator
rng('default')

% auxiliar functions
W2dBm = @(x) 10*log10(x/1e-3); % Watt to dBm
dBm2W = @(x) 1e-3*10.^(x/10);  % dBm to Watt
Pav = @(x) mean(abs(x).^2);    % Average power
AmpNorm = @(x) x/max(abs(x));  % Normalize amplitude (used for plotting)
% berqam = @(CS, SNR) 4/log2(CS)*(1 - 1/sqrt(CS))*qfunc(sqrt(3/(CS-1)*10.^(SNR/10)));
berqam = @(CS, SNR) berawgn(SNR-10*log10(log2(CS)), 'qam', CS);

% -----------------------------------------------------------------
% Define parameters
% -----------------------------------------------------------------
% Constants
c = 299792458;              % speed of light (m/s)
h = 6.626e-34;              % Planck's constant (J*s)

% Transmitter parameters (for values chosen, Butterworth filter, not interpolator, is bandwidth-limiting element)
lambdanm = 1550;            % Wavelength (nm)
Plaser = 4.4e-3;            % Average laser power (W)
kappa = 1;                  % dc-slope
interp_length = 9;          % length parameter for interpolator at Tx   
interp_cutoff = 0.8;        % cutoff parameter for interpolator at Tx
fnl = 25e9;                 % cutoff frequency of the laser
Hl = @(f) 1./(1 + 2*1j*f/fnl - (f/fnl).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
                                              % assumed 2nd-order filter
                                              % with damping = 1

% Hl = @(f) ones(size(f));


% Modulation parameters
Nc = 64;                    % total number of subcarriers = FFT size (must be even)
Nu = 52/2;                % number of used (nonzero, including complex conjugate subcarriers) subcarriers (must be even, and at most Nc/2-1)
CS = 64;                    % constellation size (4 corresponds to QPSK)
Rb = 100e9;                 % bit rate (total over two polarizations) (b/s)
Ms = Nc/Nu;                 % oversampling ratio
Rs = 2*Rb/log2(CS);           % symbol rate (Hz)

% Compute cyclic prefix length (assumes 2nd order filter with unit damping)
% Does not include CP overhead or oversampling. 
frac_incl = 0.999;                          % fraction of energy to be included within CP length
tcp = 1/Rs*(0:0.5:10);                     % vector of sampling intervals
htot = (2*pi*fnl)^2*tcp.*exp(-2*pi*fnl*tcp); % total impulse response
en_frac = cumsum(htot.^2)/sum(htot.^2);     % fraction of energy in impulse response vs. sampling interval
L = sum(en_frac<frac_incl)+1;               % number of samples that include the desired fraction of energy
Npre = L-1;                                 % minimum cyclic prefix length is L-1
                                        
CPpendB = 10*log10((Nc+Npre)/Nc);           % cyclic prefix penalty (electrical dB)

% Receiver parameters
R = 1;                      % responsivity

% % Loss and length
% alphadBkm = 0.35;       % fiber loss (dBm/km)
% conn_lossdB = 2;        % connectors losses (dB)
% Lkm = 2;                % length per fiber span (km)
% 
% % Chromatic dispersion
% Dpsnmkm = 0;            % CD (ps/(nm*km))

% Simulation parameters
Nsymb = 2^14;                   % number of OFDM symbols
Mct = 2;                        % oversampling rate for emulation of continuous time
SNR = 15:32;

% -----------------------------------------------------------------
% Compute some parameters and frequency responses
% (these calculations are done just once per simulation)
% -----------------------------------------------------------------

% Cyclic prefixes
Npre_os = ceil(Npre*Ms);        % total cyclic prefix length (oversampled) (needs to define Npre)                        
Npos_os = round(Npre_os/2);     % positive cyclic prefix length (oversampled)
Nneg_os = Npre_os - Npos_os;    % negative cyclic prefix length (oversampled)

% Time and frequency scales
fs = Rs*(Nc + Npre_os)/Nu;                                        % chip rate (Hz)
fsct = Mct*fs;                                                      % sampling frequency to emulate continuous time (Hz)
dt = 1/fsct;                                                        % time increment in emulating continuous time (s)
Ntot = Mct*(Nu*Ms + Npre_os)*Nsymb;                                 % total number of points simulated in continuous time
t = dt*(0:Ntot-1);                                                  % continuous time scale
df = fsct/Ntot;                                                     % frequency increment in continuous time (Hz)                    
f = -fsct/2:df:fsct/2-df;                                           % frequency range in continuous time                    
fc = fs/Nc*(1:2:Nu);                                                % frequency at which subcarriers are located 
tsamp = t(1:Mct:end);                                               % sampling times

% Transmitter

% Generate OFDM signal at chip rate (done in DSP)
rng('shuffle');                                         % Reinitialize the random number generator used by rand, randi, and randn with a seed based on the current time
dataTX = randi([0 CS-1], [ceil(Nu/2) Nsymb]);           % data to be encoded (symbols columnwise)
dataTXm = qammod(dataTX, CS, 0, 'gray');                % encoded QAM symbols to be modulated onto subcarriers
dataTXm = dataTXm/sqrt(2*(CS-1)/3);                     % Normalize constellation to unit energy (i.e., E(|X_n|^2) = 1)

% Perform equalization with Hl(fc)^-1
for kk = 1:Nsymb
    dataTXm(:, kk) = dataTXm(:, kk)./(kappa*Hl(fc.')); 
end

dataTXm2 = zeros(Nu, Nsymb);
dataTXm2(1:2:end, :) = dataTXm;

sig2 = 2*sum(abs(1./(kappa*Hl(fc))).^2);

% zero-pad and ensure Hermitian symmetry
% -f(N) ... -f(1) f(0) f(1) ... f(N-1) => (0 ... 0, X(-n)*, 0, X(n), 0 ... 0)
Xn = ifftshift([zeros((Nc-2*Nu)/2, Nsymb); ...     % zero-pad neg. frequencies
              flipud(conj(dataTXm2));...            % data*(-n) (Nu) 
              zeros(1, Nsymb); ...                % 0 at f == 0 (1)
              dataTXm2; ...                         % data(n)  (Nu)
              zeros((Nc-2*Nu)/2 - 1, Nsymb)], 1);   % zero-pad pos. frequencies

% Perform ifft (Nc-IDFT) columnwise
xn = Nc*ifft(Xn, Nc, 1); 

% Insert cyclic prefix
xncp = [xn(end-Npre_os+1:end, :); xn]; % insert cyclic prefix

% Parallel to serial
xt = reshape(xncp, 1, (Nc + Npre_os)*Nsymb); % time-domain ofdm signal w/ cyclic prefix

% Interpolate OFDM signal to emulate continuous-time laser drive signal (D/A conversion)
xt = interp(xt, Mct, interp_length, interp_cutoff);   

%% Clip then filter
xtc = xt;
xtc(xt < 0) = 0;
Pcf = kappa*real(ifft(ifftshift(Hl(f).*fftshift(fft(xtc)))));  % Note: real() removes residual imag 
                                                                   % part that appears due to numerical error  
                                                                 
% Detect and sample signal
Itcf = 2*R*Pcf;                      % detect signal
Itcf = Itcf(1:Mct:end);                     % Resample signal once per chip                  
Icpcf = reshape(Itcf, Nc + Npre_os, Nsymb);     % reshape into matrix form

% Remove cyclic prefix
Icf = [Icpcf(Npre_os+1:end-Nneg_os, :); Icpcf(Npos_os+1:Npre_os, :)];
%     I = circshift(Icp(Npos_os+1:end-Nneg_os, :), -Nneg_os);              % remove cyclic prefix and perform a cyclic shift
                                                                           % Note: the cyclic shift is not needed in a real receiver.
%                                                                          % It corresponds to a phase shift on each carrier, which the
%                                                                          % adaptive equalizer will compensate automatically.

% Demodulate symbols
Yncf = fft(Icf, Nc, 1)/Nc;                  % demodulated subcarrier amplitudes

%% filter then clip
Pfc = kappa*real(ifft(ifftshift(Hl(f).*fftshift(fft(xt)))));  % Note: real() removes residual imag 
                                                                   % part that appears due to numerical error 
Pfc(Pfc < 0) = 0;

% Detect and sample signal
Itfc = 2*R*Pfc;                      % detect signal
Itfc = Itfc(1:Mct:end);                     % Resample signal once per chip                  
Icpfc = reshape(Itfc, Nc + Npre_os, Nsymb);     % reshape into matrix form

% Remove cyclic prefix
Ifc = [Icpfc(Npre_os+1:end-Nneg_os, :); Icpfc(Npos_os+1:Npre_os, :)];
%     I = circshift(Icp(Npos_os+1:end-Nneg_os, :), -Nneg_os);              % remove cyclic prefix and perform a cyclic shift
                                                                           % Note: the cyclic shift is not needed in a real receiver.
%                                                                          % It corresponds to a phase shift on each carrier, which the
%                                                                          % adaptive equalizer will compensate automatically.

% Demodulate symbols
Ynfc = fft(Ifc, Nc, 1)/Nc;                  % demodulated subcarrier amplitudes


figure
stem(-Nc/2:Nc/2-1, mean(fftshift(abs(Yncf).^2), 2), '.')
hold on
stem(-Nc/2:Nc/2-1, mean(fftshift(abs(Ynfc).^2), 2), 'og')
ylabel('E(|Y_n|^2)')
xlabel('Subcarrier')
legend('Clip then filter', 'Filter then clip')
axis([-Nc/2-1 Nc/2+1 0 1.5])

