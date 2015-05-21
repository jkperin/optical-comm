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
Exc = @(sig, r) sig*(r + 1/sqrt(2*pi)*exp(-0.5*r.^2));                                 % mean value E(x+(t))
Exc2 = @(sig, r) sig^2*((r.^2 + 1).*(1 - qfunc(r)) + r/sqrt(2*pi).*exp(-r.^2/2));    % power E[x+(t)^2]

% -----------------------------------------------------------------
% Define parameters
% -----------------------------------------------------------------
% Constants
c = 299792458;              % speed of light (m/s)
h = 6.626e-34;              % Planck's constant (J*s)

% Transmitter parameters (for values chosen, Butterworth filter, not interpolator, is bandwidth-limiting element)
clip_ratio = 1;             % clip_ratio = clip_level/std, std = sqrt(P);
lambdanm = 1550;            % Wavelength (nm)
kappa = 1;                  % dc slope
interp_length = 9;          % length parameter for interpolator at Tx   
interp_cutoff = 0.8;        % cutoff parameter for interpolator at Tx
fnl = 20e9;                 % cutoff frequency of the laser
Hl = @(f) 1./(1 + 2*1j*f/fnl - (f/fnl).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
                                              % assumed 2nd-order filter
                                              % with damping = 1
                                              
Hl = @(f) ones(size(f));
% Modulation parameters
Nc = 64;                    % total number of subcarriers = FFT size (must be even)
Nu = 52;                    % number of used (nonzero, including complex conjugate subcarriers) subcarriers (must be even)
CS = 16;                    % constellation size (4 corresponds to QPSK)
Rb = 100e9;                 % bit rate (total over two polarizations) (b/s)
Ms = Nc/Nu;                 % oversampling ratio
Rs = Rb/log2(CS);           % symbol rate (Hz)

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

% Simulation parameters
Nsymb = 2^15;                   % number of OFDM symbols
Mct = 4;                        % oversampling rate for emulation of continuous time

% -----------------------------------------------------------------
% Compute some parameters and frequency responses
% (these calculations are done just once per simulation)
% -----------------------------------------------------------------

% Cyclic prefixes
Npre_os = ceil(Npre*Ms);        % total cyclic prefix length (oversampled) (needs to define Npre)                        
Npos_os = round(Npre_os/2);     % positive cyclic prefix length (oversampled)
Nneg_os = Npre_os - Npos_os;    % negative cyclic prefix length (oversampled)

% Time and frequency scales
fs = 2*Rs*(Nc + Npre_os)/Nu;                                        % chip rate (Hz)
fsct = Mct*fs;                                                      % sampling frequency to emulate continuous time (Hz)
dt = 1/fsct;                                                        % time increment in emulating continuous time (s)
Ntot = Mct*(Nu*Ms + Npre_os)*Nsymb;                                 % total number of points simulated in continuous time
t = dt*(0:Ntot-1);                                                  % continuous time scale
df = fsct/Ntot;                                                     % frequency increment in continuous time (Hz)                    
f = -fsct/2:df:fsct/2-df;                                           % frequency range in continuous time                    
fc = fs/Nc*(1:Nu/2);                                                % frequency at which subcarriers are located 
tsamp = t(1:Mct:end);                                               % sampling times

% Transmitter
lambda = lambdanm*1e-9;                 % optical wavelength (m)
nu = c/lambda;                          % optical frequency (Hz)

%% Generate OFDM signal at chip rate (done in DSP)
rng('shuffle');                                         % Reinitialize the random number generator used by rand, randi, and randn with a seed based on the current time
dataTX = randi([0 CS-1], [Nu/2 Nsymb]);                   % data to be encoded (symbols columnwise)
dataTXm = qammod(dataTX, CS, 0, 'gray');             % encoded QAM symbols to be modulated onto subcarriers

% Perform equalization with Hl(fc)^-1
for kk = 1:Nsymb
    dataTXm(:, kk) = dataTXm(:, kk)./(kappa*Hl(fc.'));
end

% zero-pad and ensure Hermitian symmetry
% -f(N) ... -f(1) f(0) f(1) ... f(N-1) => (0 ... 0, X(-n)*, 0, X(n), 0 ... 0)
Xn = ifftshift([zeros((Nc-Nu)/2, Nsymb); ...     % zero-pad neg. frequencies
          flipud(conj(dataTXm));...            % data*(-n) (Nu) 
          zeros(1, Nsymb); ...                % 0 at f == 0 (1)
          dataTXm; ...                         % data(n)  (Nu)
          zeros((Nc-Nu)/2 - 1, Nsymb)], 1);   % zero-pad pos. frequencies

% Perform ifft (Nc-IDFT) columnwise
xn = Nc*ifft(Xn, Nc, 1); 

% Insert cyclic prefix
xncp = [xn(end-Npre_os+1:end, :); xn]; % insert cyclic prefix

% Parallel to serial
xt = reshape(xncp, 1, (Nc + Npre_os)*Nsymb); % time-domain ofdm signal w/ cyclic prefix

% Interpolate OFDM signal to emulate continuous-time laser drive signal (D/A conversion)
xt = interp(xt, Mct, interp_length, interp_cutoff);   

% Hard clip at negative level
sig = sqrt(2*sum(abs(1./(kappa*Hl(fc))).^2))*sqrt(2*(CS-1)/3);
clip_level = clip_ratio*sig;
clip_prob = sum(xt < -clip_level)/numel(xt);
xtc = xt + clip_level;                 % add dc-bias to make signal positive
xtc(xtc < 0) = 0;

varxt = var(xt);
Extcs = mean(xtc)
Extc = Exc(sig, clip_ratio)
Extc2s = mean(abs(xtc).^2)
Extc2 = Exc2(sig, clip_ratio)

Extc = sig*(clip_ratio + 1/sqrt(2*pi)*exp(-clip_ratio^2/2));

% Apply frequency response of the laser
P = kappa*real(ifft(ifftshift(Hl(f).*fftshift(fft(xtc)))));  % Note: real() is use to remove residual imag
                                                               % part that appears due to numerical error          
% Detect and sample signal
It = R*(P - Extc*kappa*Hl(0));                      % detect signal and remove DC bias
It = It(1:Mct:end);                     % Resample signal once per chip                  

Icp = reshape(It, Nc + Npre_os, Nsymb);     % reshape into matrix form

% Remove cyclic prefix
I = [Icp(Npre_os+1:end-Nneg_os, :); Icp(Npos_os+1:Npre_os, :)];
%     I = circshift(Icp(Npos_os+1:end-Nneg_os, :), -Nneg_os);              % remove cyclic prefix and perform a cyclic shift
                                                                       % Note: the cyclic shift is not needed in a real receiver.
                                                                       % It corresponds to a phase shift on each carrier, which the
                                                                       % adaptive equalizer will compensate automatically.

% Demodulate symbols
Yn = fft(I, Nc, 1)/Nc;                  % demodulated subcarrier amplitudes
Yn(1,:) = 0;
dataRXm = Yn(2:Nu/2+1, :);            % used subcarrier amplitudes (complex conjugate subcarriers are ignored)

dataRX = qamdemod(dataRXm, CS, 0, 'gray');      % detected data    

[num_errs, ber] = biterr(reshape(dataTX, 1, Nu/2*Nsymb), reshape(dataRX, 1, Nu/2*Nsymb));   % number of bit errors and bit-error probability  


%% Figures
figure
subplot(221)
plot(sqrt(2*(CS-1)/3)*dataTXm, '.')
axis square
axis(sqrt(CS)*[-1 1 -1 1])
title('Scatterplot at Transmitter');

subplot(222)
plot(dataRXm, '.')
axis(sqrt(CS)*[-1 1 -1 1])
axis square
title('Scatterplot after clipping');

subplot(223)
stem(-Nc/2:Nc/2-1, mean(fftshift(abs(Yn).^2), 2), 'filled')
ylabel('E(|Y_n|^2)')
xlabel('Subcarrier')
% axis([-Nc/2-1 Nc/2+1 0 1.5])

subplot(224)
axis off
text(-0.1,1, sprintf('R_{s} = %.1f GHz, f_{nl} = %.1f GHz', Rs/1e9, fnl/1e9));
text(-0.1,0.8, sprintf('N = %d, CS = %d, Nu = %d', Nc, CS, Nu));
text(-0.1,0.6, sprintf('L = %d, CP Penalty = %.2f dB', L, CPpendB));
text(-0.1,0.4, sprintf('Clip. Ratio = %d, P_{clip} = %.2G (%.2G Gaussian)', clip_ratio, mean(clip_prob), qfunc(clip_ratio)));
% text(-0.1,0.2, sprintf('Estimated penalty = %.2f dB @ BER 10^{-12}', SNRpenalty));
% text(-0.1,0.0,[ 'f_{VCSEL} = ' num2str(fl/1e9,3) ' GHz, f_{MMF} = ' num2str(ff/1e9,3) ' GHz, Loss = ' num2str(lossdB,3) ' dB']);

% 
dtx = reshape(dataTXm, numel(dataTXm), 1);  
dtx2 = [real(dtx) imag(dtx)];

drx = reshape(dataRXm, numel(dataRXm), 1);
drx2 = [real(drx) imag(drx)];
figure
hist3(drx2, [75 75])
set(gcf,'renderer','opengl');
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');

figure
hist3(drx2, [75 75])
set(gcf,'renderer','opengl');
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
view(2)

dclip = drx2-dtx2;
figure
hist3(dclip, [50 50])
set(gcf,'renderer','opengl');
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
view(2)
axis square

sig2 = 2*(CS-1)/3;
EXY = sig2*(1-qfunc(clip_ratio));

u = dclip(:,1) + 1j*dclip(:,2);
EXu = -sig2*qfunc(clip_ratio);

var(u);
2*sig2*qfunc(clip_ratio);

r = clip_ratio;
var(drx/sqrt(2*(CS-1)/3))
fr = Exc2(1, r) - Exc(1, r)^2

% x = qammod(0:CS-1, CS);
% for k = 1:length(x)
%     u = drx(dtx == x(k)) - dtx(dtx == x(k));
%     vark(k) = var(u);
% end

clip_prob
qfunc(clip_ratio)


% dc = dclip(:,1) + 1j*dclip(:,2);


