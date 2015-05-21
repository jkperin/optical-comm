% ofdm_twopol_adaptive
% Joseph M. Kahn, September 4, 2008

clear all;
close all;

tic

% -----------------------------------------------------------------
% Define parameters
% -----------------------------------------------------------------

% Constants
j = sqrt(-1);
c = 299792458;          % speed of light (m/s)
h = 6.626e-34;          % Planck's constant (J*s)

% Modulation parameters
Nu = 52;                % number of used subcarriers (must be even)
Nc = 64;                % total number of subcarriers = FFT size (must be even)
Npre = 10;              % total cyclic prefix length (ref. to 1 sample/chip)
Ms = Nc/Nu;             % oversampling ratio
CS = 4;                 % constellation size (4 corresponds to QPSK)
Rb = 112e9;             % bit rate (total over two polarizations) (b/s)
Rs = Rb/(2*log2(CS));       % symbol rate (Hz)

% Transmitter parameters (for values chosen, Butterworth filter, not interpolator, is bandwidth-limiting element)
clip_level = 0.60;          % clipping level (V) (should adjust to get best results)
Vpi = 5;                    % switching voltage of modulator
txquantrange = 2*Vpi;       % range between max and min levels of Tx D/A (V) (should be 2*Vpi when clipping level set properly)
txquantbits = 8;            % bits of quantization of Tx D/A
interp_length = 9;          % length parameter for interpolator at Tx   
interp_cutoff = 0.8;        % cutoff parameter for interpolator at Tx
tx_but_ord = 5;             % order of Butterworth filter at Tx
tx_but_cutoff_norm = 0.5;   % cutoff frequency of Butterworth filter Tx (normalized to chip rate)
lambdanm = 1550;            % Wavelength (nm)
PlaserdBm = 6;              % power of CW laser
PtxdBm_ber = [-9:1.5:0];   % transmit signal powers at which to measure BER
PtxdBm_adapt = -3;         % power at which to do adaptation
PtxdBm_display = -3;         % power at which to display received spectrum, scatter plots of signal constellations, etc.
PtxdBm_scan = [PtxdBm_adapt, PtxdBm_ber, PtxdBm_display];

% Receiver parameters
rx_but_ord = 5;             % order of Butterworth filter at Rx
rx_but_cutoff_norm = 1/2;   % cutoff frequency of Butterworth filter Rx (normalized to chip rate)
rxagcrms = 0.15;            % filtered signal is scaled to this rms voltage before sampling and quantization
rxquantrange = 0.75;        % range between max and min levels of Rx A/D (V)
rxquantbits = 8;            % bits of quantization of Rx A/D

% LMS algorithm
muinit = 0.005;              % initial value of step size parameter
mufactor = 0.1;            % a factor less than one
muinterval = 768;            % multiply mu by mufactor each time k increments by muinterval
mushifts = 1;               % Number of gear shifts. Stop multiplying mu by mufactor once k reaches muinterval*(mushifts+1)
startdd = 128;       % symbol at which to change from trained mode to decision-directed mode

% Loss and length
alphadBkm = 0.25;       % fiber loss (dBm/km)
Nspan = 25;             % number of fiber spans
Lspankm = 80;           % length per fiber span (km)
Lkm = Nspan*Lspankm;    % fiber length (km)

% Chromatic dispersion
Dpsnmkm = 17;           % CD (ps/(nm*km))
comppercent = 96;       % CD optical compensation (%)

% Polarization-mode dispersion
pmd = 'y';
taumeanps = 30;                 % mean DGD (ps)
Nsect = 10;                     % number of sections used to generate PMD
realiznumb = 1;                 % realization number
jmoldnew = 'old';               % if 'old': use existing Jones matrix jmfilename.mat
                                % if 'new': create and store new Jones matrix jmfilename.mat
                                % Note: must choose 'new' if change oversampling rate, number of carriers or cyclic prefix length.
% EDFA
NFdB = 5.0;                     % EDFA noise figure (dB)

% ROADM
ROADMord = 2;                   % ROADM super-Gaussian order
ROADMbw = 50e9;                 % ROADM bandwidth at -3 dB (Hz)
pre_ROADM_ASE_percent = 50;     % percent of ASE addded before ROADM (unitless, between 0 and 100)

% Simulation parameters
newsignal = 'n';                % 'y': generate new signal for each power level, 'n': generate one signal for all power levels
newnoise = 'n';                 % 'n': generate new noise for each power level, 'n': generate one noise for all power levels
Nsymb = 2^10;                   % number of OFDM symbols
Mct = 2;                        % oversampling rate for emulation of continuous time
fplotmax = 25e9;                % maximum frequency for plots of frequency-dependent quantities
osnrnoisebw = 50e9;             % optical filter bandwidth for measurement of optical signal-to-noise ratio
Npdg = 8;                       % block length for periodogram
ch_sel = [1 Nu/2 Nu/2+1 Nu];    % channels selected to display scatter plot of signal constellation
smoothlength = 8;               % window length for smooothing error plots
eqmagmindB = -20;               % min and max magnitudes and phases for plots of equalizer coefficients
eqmagmaxdB = 10;
eqphsmin = -1*pi;
eqphsmax = 5*pi;

% -----------------------------------------------------------------
% Compute some parameters and frequency responses
% (these calculations are done just once per simulation)
% -----------------------------------------------------------------

% Cyclic prefixes
Npre_os = ceil(Npre*Ms);        % total cyclic prefix length (oversampled)                            
Npos_os = round(Npre_os/2);     % positive cyclic prefix length (oversampled)
Nneg_os = Npre_os - Npos_os;    % negative cyclic prefix length (oversampled)

% Time and frequency scales
fs = Rs*(Nc + Npre_os)/Nu;                                          % chip rate (Hz)
fsct = Mct*fs;                                                      % sampling frequency to emulate continuous time (Hz)
dt = 1/fsct;                                                        % time increment in emulating continuous time (Hz)
Ntot = Mct*(Nu*Ms + Npre_os)*Nsymb;                                 % total number of points simulated in continuous time
t = dt*(0:Ntot-1);                                                  % continuous time scale
tsamp = t(1:Mct:end);                                               % sampling times
deltaomega = 2*pi*fsct/Ntot;                                        % increment of continuous frequency
omega = deltaomega*[0:(Ntot/2)-1 -Ntot/2:-1];                       % continuous frequency scale
plotind = find(abs(omega)/2/pi<=fplotmax);                          % indices of continuous frequencies to plot
[omegaplotsort, plotindsort] = sort(omega(plotind));
r = Ntot/(Nc*Mct);
ucind = [1/2*r+1:r:(Nu-1)/2*r+1 (-Nu+1)/2*r+Ntot+1:r:-1/2*r+Ntot+1];% indices of used subcarriers in continuous frequency scale
[omegaucsort, ucindsort] = sort(omega(ucind));

% Carrier to center spectrum at Tx and uncenter it at Rx
vcent = exp(+j*2*pi*fs/(2*Nc)*(0:1/fs:((Nc+Npre_os)*Nsymb-1)/fs));        % modulated carrier to center spectrum

% Butterworth filter for Tx
[num,denom] = butter(tx_but_ord,tx_but_cutoff_norm*fs*2*pi*(Nu+1)/Nc,'s'); 
Htxlpf = polyval(num, j*omega)./polyval(denom, j*omega);
gdtxlpf = -diff(unwrap(angle(Htxlpf)))/deltaomega;              
gdtx = gdtxlpf(1);
Htxlpf = Htxlpf.*exp(j*omega*gdtx);                             % continuous frequency response of Tx filter (delay removed)

% Butterworth filter for Rx
[num,denom] = butter(rx_but_ord,rx_but_cutoff_norm*fs*2*pi*(Nu+1)/Nc,'s');
Hrxlpf = polyval(num, j*omega)./polyval(denom, j*omega);
gdrxlpf = -diff(unwrap(angle(Hrxlpf)))/deltaomega;
gdrx = gdrxlpf(1);
Hrxlpf = Hrxlpf.*exp(j*omega*gdrx);                             % continuous frequency response of Rx filter (delay removed)

% Transmitter
lambda = lambdanm*1e-9;                 % optical wavelength (m)
nu = c/lambda;                          % optical frequency (Hz)
Plaser = 1e-3*10^(PlaserdBm/10);        % laser power (W)
Ptx_scan = 1e-3*10.^(PtxdBm_scan/10);    % transmitted power (W)

% Chromatic dispersion
Dsmm = 1e-6*Dpsnmkm;                            % convert units of D
beta2 = -lambda^2/(2*pi*c)*Dsmm;                % compute beta2
Lm = 1e3*Lkm;                                   % convert units of L
Luncomp = Lm*(100-comppercent)/100;             % equivalent uncompensated length of fiber
Hcd = exp(-j*1/2*Luncomp*beta2*omega.^2);       % continous frequency response of fiber CD

% Polarization-mode dispersion
if pmd == 'y';
    jmfilename = ['jm' sprintf('%d',taumeanps) 'ps'...
        sprintf('%d',Nsect) 'sect' sprintf('%d',realiznumb)]; % file name for Jones matrix

    if all(jmoldnew == 'old')

        eval(['load ' jmfilename]);

    elseif all(jmoldnew == 'new')

        dtau = taumeanps*1e-12/sqrt(Nsect);

        phi = rand(1, 3)*2*pi;

        U1 = [exp(-j*phi(1)/2), 0; 0 exp(j*phi(1)/2)];
        U2 = [cos(phi(2)/2) -j*sin(phi(2)/2); -j*sin(phi(2)/2) cos(phi(2)/2)];
        U3 = [cos(phi(3)/2) -sin(phi(3)/2); sin(phi(3)/2) cos(phi(3)/2)];

        U = U1*U2*U3;

        M = zeros(2,2,length(omega));
        for m = 1:length(omega);
            M(:,:,m) = U;
        end

        for k = 1:Nsect
            phi = rand(1, 3)*2*pi;

            U1 = [exp(-j*phi(1)/2), 0; 0 exp(j*phi(1)/2)];
            U2 = [cos(phi(2)/2) -j*sin(phi(2)/2); -j*sin(phi(2)/2) cos(phi(2)/2)];
            U3 = [cos(phi(3)/2) -sin(phi(3)/2); sin(phi(3)/2) cos(phi(3)/2)];

            U = U1*U2*U3;

            for m = 1:length(omega)
                D = [exp(-j*dtau*omega(m)/2), 0; 0, exp(j*dtau*omega(m)/2)];
                M(:,:,m) = M(:,:,m)*U'*D*U;
            end
        end

        eval(['save ' jmfilename ' M']);
    end
    
elseif pmd == 'n'

M = zeros(2,2,length(omega));
M(1,1,:) = ones(size(omega));
M(2,2,:) = ones(size(omega));
        
end

% ROADM
b = log(2)/2*(pi*ROADMbw)^(-2*ROADMord);    % parameter describing super-Gaussian freq. resp.
Hroadm = exp(-b*omega.^(2*ROADMord));       % super-Gaussian freq. resp.

% Optical amplification
NF = 10^(NFdB/10);                  % noise figure
nsp = NF/2;                         % inversion factor
G = 10^(alphadBkm*Lspankm/10);      % gain per span, assumed to equal fiber loss per span
Ssp = nsp*(G-1)*h*nu;               % PSD of ASE from one amplifier (one-sided PSD of real passband noise in one pol.,
                                    % also equal to two-sided PSD of complex baseband noise in one pol.)
Spreroadm = pre_ROADM_ASE_percent/100*Nspan*Ssp;        % PSD of the total ASE added before ROADM
Spostroadm = (1-pre_ROADM_ASE_percent/100)*Nspan*Ssp; 	% PSD of the total ASE added after ROADM

% OSNR
% This is computed including two signal polarizations and two noise polarizations over an optical bandwidth
% osnrnoisebw. This does not take into account that noise introduced before ROADM is filtered. If the ROADM
% bandwidth is less than the OSNR measurement bandwidth, the actual OSNR will be larger than that estimated here. 
OSNR_scan = Ptx_scan/(2*Nspan*Ssp*osnrnoisebw);
OSNR_dB_scan = 10*log10(OSNR_scan);
OSNR_dB_adapt = OSNR_dB_scan(1);
OSNR_dB_display = OSNR_dB_scan(end);
OSNR_dB_ber = OSNR_dB_scan(2:end-1);

% Ideal (Non-Adaptive) Equalizer
Hs = 1./(Hcd(ucind).*Hroadm(ucind).*Hrxlpf(ucind));        % frequency response of scalar equalizer (for CD/ROADM/Rx filter) at used subcarrier frequencies  
Heq = zeros(2,2,Nu);
for m = 1:Nu;
    Heq(:,:,m) = Hs(m)*inv(M(:,:,ucind(m)));                      % frequency response of matrix equalizer (for CD/ROADM/Rx filter/PMD) at used subcarrier frequencies
end

% LMS Algorithm
muvec = muinit*mufactor.^( min((mushifts+1),ceil((1:Nsymb)/muinterval))-1);     % vector of step size parameters
decwtvec = (1:Nsymb)<startdd;   % vector of decision weights for LMS algorithm (1 when use training symbols, 0 when use decision symbols)

% -----------------------------------------------------------------
% Perform simulation of signal waveform and its propagation
% (these calculations are done at several different optical power levels)
% -----------------------------------------------------------------

% Note: "1" and "2" denote x and y polarizations, respectively

% Initialize some vectors that are recorded at different powers
num_errs_scan = [];                                     % vector of number of bit errors for fixed equalizer
ber_scan = [];                                          % vector of bit-error probabilities for fixed equalizer
num_errs_scan_ad = [];                                     % vector of number of bit errors for adaptive equalizer
ber_scan_ad = [];                                          % vector of bit-error probabilities for adaptive equalizer

% Scan the power
for k = 1:length(PtxdBm_scan);
    
    if (newsignal == 'y') || (k == 1);
        
        % Generate OFDM signal at chip rate (done in DSP)
        d1 = randint(Nu,Nsymb,CS,sum(100*clock));                           % data to be encoded
        d2 = randint(Nu,Nsymb,CS,sum(100*clock));                           % data to be encoded
        V1 = qammod(d1,CS,0,'gray');                                        % encoded QAM symbols to be modulated onto subcarriers
        V2 = qammod(d2,CS,0,'gray');                                        % encoded QAM symbols to be modulated onto subcarriers
        Vscale1 = V1.*(1./Htxlpf(ucind).'*ones(1,Nsymb));                   % scale the subcarriers to compensate for Tx filter
        Vscale2 = V2.*(1./Htxlpf(ucind).'*ones(1,Nsymb));                   % scale the subcarriers to compensate for Tx filter
        Vz1 = [Vscale1(1:Nu/2,:); zeros((Nc-Nu),Nsymb);  Vscale1(Nu/2+1:Nu,:)];             % insert zero subcarriers to implement oversampling by Ms = Nc/Nu
        Vz2 = [Vscale2(1:Nu/2,:); zeros((Nc-Nu),Nsymb);  Vscale2(Nu/2+1:Nu,:)];             % insert zero subcarriers to implement oversampling by Ms = Nc/Nu
        vz1 = fftshift(ifft(Vz1,[],1),1);                                                   % time-domain waveform
        vz2 = fftshift(ifft(Vz2,[],1),1);                                                   % time-domain waveformvcp = [vz(end-Npre_os+1:1:end,:); vz];                                          % time-domain waveform with cyclic prefix samples prepended
        vcp1 = [vz1(end-Npre_os+1:end,:); vz1];                                             % time-domain waveform with cyclic prefix samples prepended
        vcp2 = [vz2(end-Npre_os+1:end,:); vz2];                                             % time-domain waveform with cyclic prefix samples prepended
        vcpser1 = reshape(vcp1,1,(Nc+Npre_os)*Nsymb);                                       % time-domain waveform serialized
        vcpser2 = reshape(vcp2,1,(Nc+Npre_os)*Nsymb);                                       % time-domain waveform serialized
        vcpcent1 = vcpser1.*vcent;                                                          % time-domain waveform with spectrum centered
        vcpcent2 = vcpser2.*vcent;                                                          % time-domain waveform with spectrum centered

        % Perform nonlinear predistortion of waveform to compensate for MZ modulator (done in DSP)
        vcpcent1i = real(vcpcent1);                             % in-phase component
        vcpcent1q = imag(vcpcent1);                             % quadrature component
        vcpcent2i = real(vcpcent2);                             % in-phase component
        vcpcent2q = imag(vcpcent2);                             % quadrature component
        vclip1i = max(-clip_level,min(clip_level,vcpcent1i));   % in-phase after clipping
        vclip1q = max(-clip_level,min(clip_level,vcpcent1q));   % quadrature after clipping
        vclip2i = max(-clip_level,min(clip_level,vcpcent2i));   % in-phase after clipping
        vclip2q = max(-clip_level,min(clip_level,vcpcent2q));   % quadrature after clipping
        vdist1i = 2*Vpi/pi*asin(vclip1i/clip_level);            % in-phase after predistortion
        vdist1q = 2*Vpi/pi*asin(vclip1q/clip_level);            % quadrature after predistortion
        vdist2i = 2*Vpi/pi*asin(vclip2i/clip_level);            % in-phase after predistortion
        vdist2q = 2*Vpi/pi*asin(vclip2q/clip_level);            % quadrature after predistortion
        vquant1i = quantz(vdist1i,txquantrange,txquantbits);    % in-phase including effect of quantization
        vquant1q = quantz(vdist1q,txquantrange,txquantbits);    % quadrature including effect of quantization
        vquant2i = quantz(vdist2i,txquantrange,txquantbits);    % in-phase including effect of quantization
        vquant2q = quantz(vdist2q,txquantrange,txquantbits);    % quadrature including effect of quantization

        % Interpolate OFDM signal to emulate continuous-time modulator drive signal (D/A conversion)
        vinterp1i = interp(vquant1i,Mct,interp_length,interp_cutoff);   % in-phase after nearly ideal interpolation
        vinterp1q = interp(vquant1q,Mct,interp_length,interp_cutoff);   % quadrature after nearly ideal interpolation
        vinterp2i = interp(vquant2i,Mct,interp_length,interp_cutoff);   % in-phase after nearly ideal interpolation
        vinterp2q = interp(vquant2q,Mct,interp_length,interp_cutoff);   % quadrature after nearly ideal interpolation
        vmod1i = real(ifft(fft(vinterp1i).*Htxlpf));                    % in-phase after Butterworth filter
        vmod1q = real(ifft(fft(vinterp1q).*Htxlpf));                    % quadrature after Butterworth filter
        vmod2i = real(ifft(fft(vinterp2i).*Htxlpf));                    % in-phase after Butterworth filter
        vmod2q = real(ifft(fft(vinterp2q).*Htxlpf));                    % quadrature after Butterworth filter

        % Perform modulation using MZ modulator
        emod1i = sqrt(Plaser)/2*sin(pi*vmod1i/(2*Vpi)); % in-phase modulator output
        emod1q = sqrt(Plaser)/2*sin(pi*vmod1q/(2*Vpi)); % quadrature modulator output
        emod2i = sqrt(Plaser)/2*sin(pi*vmod2i/(2*Vpi)); % in-phase modulator output
        emod2q = sqrt(Plaser)/2*sin(pi*vmod2q/(2*Vpi)); % quadrature modulator output
        emod1 = emod1i + j*emod1q;                      % modulator output
        emod2 = emod2i + j*emod2q;                      % modulator output
        Pmod = mean(abs(emod1).^2 + abs(emod2).^2);     % average power of modulator output
        eta = Pmod/Plaser;                              % modulation efficiency (W/W)
        etadB = 10*log10(eta);                          % modulation efficiency (dB)
   
    end
    
    % Optically amplify transmitted signal to desired power level
    Ptx = Ptx_scan(k);      % Transmitted power (W)
    Gtx = Ptx/Pmod;         % Amplifier gain required to reach desired transmit power level (W/W)
    GtxdB = 10*log10(Gtx);  % Amplifier gain (dB)
    etx1 = sqrt(Gtx)*emod1;   % transmitted optical signal
    etx2 = sqrt(Gtx)*emod2;   % transmitted optical signal

    % Chromatic dispersion
    ecd1 = ifft(fft(etx1).*Hcd);
    ecd2 = ifft(fft(etx2).*Hcd);
    
    % Polarization-mode dispersion
    epmd1 = ifft(fft(ecd1).*reshape(M(1,1,:),size(omega)) + fft(ecd2).*reshape(M(1,2,:),size(omega)));
    epmd2 = ifft(fft(ecd1).*reshape(M(2,1,:),size(omega)) + fft(ecd2).*reshape(M(2,2,:),size(omega)));
    
    if (newnoise == 'y') || (k == 1);
    
        % Generate ASE before ROADM
        nb1i = sqrt(Spreroadm/2/dt)*randn(1,Ntot);      % I noise with variance Spreroadm/2/dt
        nb1q = sqrt(Spreroadm/2/dt)*randn(1,Ntot);      % Q noise with variance Spreroadm/2/dt
        nb2i = sqrt(Spreroadm/2/dt)*randn(1,Ntot);      % I noise with variance Spreroadm/2/dt
        nb2q = sqrt(Spreroadm/2/dt)*randn(1,Ntot);      % Q noise with variance Spreroadm/2/dt
        nb1 = nb1i + j*nb1q;                            % complex baseband noise with variance Spreroadm/dt
        nb2 = nb2i + j*nb2q;                            % complex baseband noise with variance Spreroadm/dt
    
    end

    % Optical amplification before ROADM
    eampb1 = epmd1 + nb1;
    eampb2 = epmd2 + nb2;                          

    % ROADMs
    eroadm1 = ifft(fft(eampb1).*Hroadm);
    eroadm2 = ifft(fft(eampb2).*Hroadm);    
    
    if (newnoise == 'y') || (k == 1);

        % Generate ASE after ROADM
        na1i = sqrt(Spostroadm/2/dt)*randn(1,Ntot);     % I noise with variance Spostroadm/2/dt
        na1q = sqrt(Spostroadm/2/dt)*randn(1,Ntot);     % Q noise with variance Spostroadm/2/dt
        na2i = sqrt(Spostroadm/2/dt)*randn(1,Ntot);     % I noise with variance Spostroadm/2/dt
        na2q = sqrt(Spostroadm/2/dt)*randn(1,Ntot);     % Q noise with variance Spostroadm/2/dt
        na1 = na1i + j*na1q;                            % complex baseband noise with variance Spostroadm/dt
        na2 = na2i + j*na2q;                            % complex baseband noise with variance Spostroadm/dt
    
    elseif newnoise == 'n'
    end
    
    % Optical amplification after ROADM
    eampa1 = eroadm1 + na1;
    eampa2 = eroadm2 + na2;    
    
    % Final reception
    erx1 = eampa1;
    erx2 = eampa2;
    
    % Homodyne downconversion
    vrx1 = erx1;          % downconverted receiver output
    vrx2 = erx2;          % downconverted receiver output

    % Filter, AGC, sample, quantize
    vfilt1 = ifft(fft(vrx1).*Hrxlpf);                                               % output of anti-aliasing lowpass filter
    vfilt2 = ifft(fft(vrx2).*Hrxlpf);                                               % output of anti-aliasing lowpass filter
    Aagc = rxagcrms/(0.5*(sqrt(mean(abs(vfilt1).^2))+sqrt(mean(abs(vfilt2).^2))));  % gain of AGC
    vagc1 = Aagc*vfilt1;                        % output of AGC
    vagc2 = Aagc*vfilt2;                        % output of AGC
    vsamp1 = vagc1(1:Mct:end);                                                    	% sample once per chip
    vsamp2 = vagc2(1:Mct:end);                                                    	% sample once per chip
    vquant1 = quantz(real(vsamp1),rxquantrange,rxquantbits) + ...
        j*quantz(imag(vsamp1),rxquantrange,rxquantbits);                            % quantized samples 
    vquant2 = quantz(real(vsamp2),rxquantrange,rxquantbits) + ...
        j*quantz(imag(vsamp2),rxquantrange,rxquantbits);                            % quantized samples 
    
    % Undo spectral centering, remove cyclic prefix
    vfiltuncent1 = vquant1.*conj(vcent);                                            % spectral centering undone
    vfiltuncent2 = vquant2.*conj(vcent);                                            % spectral centering undone
    vreshape1 = reshape(vfiltuncent1,Nc+Npre_os,Nsymb);                         	% reshape into matrix
    vreshape2 = reshape(vfiltuncent2,Nc+Npre_os,Nsymb);                           	% reshape into matrix
    vnoprefix1 = circshift(vreshape1(Npos_os+1:end-Nneg_os,:),Nc/2-Nneg_os);      	% remove cyclic prefix and perform a cyclic shift
    vnoprefix2 = circshift(vreshape2(Npos_os+1:end-Nneg_os,:),Nc/2-Nneg_os);       	% remove cyclic prefix and perform a cyclic shift
                                                                                    % Note: the cyclic shift is not needed in a real receiver.
                                                                                    % It corresponds to a phase shift on each carrier, which the
                                                                                    % adaptive equalizer will compensate automatically.

    % Demodulate symbols
    y1 = fft(vnoprefix1,[],1);                  % demodulated subcarrier amplitudes
    y2 = fft(vnoprefix2,[],1);                  % demodulated subcarrier amplitudes
    yu1 = y1([1:Nu/2,Nc-Nu/2+1:Nc],:);          % used subcarrier amplitudes
    yu2 = y2([1:Nu/2,Nc-Nu/2+1:Nc],:);          % used subcarrier amplitudes
    
    % Adapt the adaptive equalizer
    if k == 1;
        
        % Initialize equalizer and error
        Heq_ad = zeros(2,2,Nu);                 % initialize equalizer  
        epsilon = zeros(2,Nu,Nsymb);            % initialize error (based on transmitted symbols)
        
        % Adapt equalizer
         for l = 1:Nsymb;            
             for m = 1:Nu;
                Zad_m_l = Heq_ad(:,:,m)*[yu1(m,l);yu2(m,l)];                                          % equalizer output
                zad_m_l = qammod(qamdemod(Zad_m_l,CS,0,'gray'),CS,0,'gray');                           % detected data
                epsilon_D_m_l = Zad_m_l - zad_m_l;                                                      % error based on decision symbols
                epsilon_m_l = Zad_m_l - [V1(m,l);V2(m,l)];                                              % error based on transmitted symbols
                epsilon(:,m,l) = epsilon_m_l;                                                           % record error based on transmitted symbols for plotting
                epsilon_used_m_l = decwtvec(l)*epsilon_m_l + ~decwtvec(l)*epsilon_D_m_l;                % error used in LMS algorithm
                Heq_ad(:,:,m) = Heq_ad(:,:,m) - 2*muvec(l)*kron(epsilon_used_m_l,[yu1(m,l);yu2(m,l)]');	% LMS update equation
            end
         end
    end
     
    % Perform equalization, detection, error counting
    if ~(k==1);
        
        % Perform equalization
        % fixed equalizer
        Z1 = zeros(Nu,Nsymb);
        Z2 = zeros(Nu,Nsymb);
        for m = 1:Nu;
            Z1(m,:) = yu1(m,:)*Heq(1,1,m) + yu2(m,:)*Heq(1,2,m);   % equalized symbol sequences on used subcarriers
            Z2(m,:) = yu1(m,:)*Heq(2,1,m) + yu2(m,:)*Heq(2,2,m);   % equalized symbol sequences on used subcarriers
        end  
        % adaptive equalizer
        Z1ad = zeros(Nu,Nsymb);
        Z2ad = zeros(Nu,Nsymb);
        for m = 1:Nu;
            Z1ad(m,:) = yu1(m,:)*Heq_ad(1,1,m)+yu2(m,:)*Heq_ad(1,2,m);   % equalized symbol sequences on used subcarriers
            Z2ad(m,:) = yu1(m,:)*Heq_ad(2,1,m)+yu2(m,:)*Heq_ad(2,2,m);   % equalized symbol sequences on used subcarriers
        end
                
        % perform detection
        % fixed equalizer     
        z1 = qamdemod(Z1,CS,0,'gray');                 % detected data
        z2 = qamdemod(Z2,CS,0,'gray');                 % detected data
        % adaptive equalizer
        z1ad = qamdemod(Z1ad,CS,0,'gray');          	% detected data
        z2ad = qamdemod(Z2ad,CS,0,'gray');          	% detected data
        
        % Count errors
        if ~((k==1)||(k==length(PtxdBm_scan)))
            % fixed equalizer
            [num_errs,ber] = biterr(reshape([d1;d2],1,2*Nu*Nsymb),reshape([z1;z2],1,2*Nu*Nsymb));   % number of bit errors and bit-error probability
            num_errs_scan = [num_errs_scan, num_errs];
            ber_scan = [ber_scan, ber];
            % adaptive equalizer
            [num_errs_ad,ber_ad] = biterr(reshape([d1;d2],1,2*Nu*Nsymb),reshape([z1ad;z2ad],1,2*Nu*Nsymb));   % number of bit errors and bit-error probability
            num_errs_scan_ad = [num_errs_scan_ad, num_errs_ad];
            ber_scan_ad = [ber_scan_ad, ber_ad];
        end
        
    end   
    
end

% -----------------------------------------------------------------
% Plot PMD, spectra, waveforms and samples, constellations, error probability
% -----------------------------------------------------------------

% % Polarization mode dispersion
% if (pmd == 'y')
%     % Differential group delay
%     for m = 1:length(omega)-1;
%         tau(m) = 2/(omega(2)-omega(1))*sqrt(det(M(:,:,m+1)-M(:,:,m)));
%     end
%     tauplot = tau(plotind(1:length(plotind)-1));
%     tauplotsort = tauplot(plotindsort(find(plotindsort<=length(plotind)-1)));
%     figure
%     clf
%     subplot(211)
%     plot(omegaplotsort(1:length(plotind)-1)/2/pi/1e9, abs(tauplotsort)/1e-12, 'b-');
%     h=get(gca,'children'); set(h,'linewidth',1.5)
%     xlabel('Frequency \omega/2\pi (GHz)');
%     ylabel('|\Delta\tau(\omega)| (ps)');
%     title([...
%         %     ['\Delta\tau_{mean} = ' num2str(taumeanps,3) ' ps'] ', '...
%         %     ['\itN\rm_{sections} = ' num2str(Nsect,3)] ', '...
%         %     ['realization number = ' num2str(realiznumb)] ': '...
%         ['Differential group delay vs. \omega']]);
% 
%     % Output SOP
%     s = zeros(length(plotind), 3);
%     for m = plotind
%         sx = M(1,1,m);
%         sy = M(2,1,m);
%         s1(m) = sx*conj(sx) - sy*conj(sy);
%         s2(m) = sx*conj(sy) + conj(sx)*sy;
%         s3(m) = j*(sx*conj(sy) - conj(sx)*sy);
%     end
%     subplot(212)
%     plot3(s1, s2, s3,'b-');
%     hold on
%     plot3(1, 0, 0, 'r*');
%     h=get(gca,'children'); set(h,'linewidth',1.5)
%     [Xs,Ys,Zs] = sphere(100);
%     surf(Xs,Ys,Zs);
%     shading interp;
%     colormap('gray')
%     alpha(1)
%     axis image
%     view(135,30);
%     rotate3d on;
%     xlabel('S_{1}/S_{0}')
%     ylabel('S_{2}/S_{0}')
%     zlabel('S_{3}/S_{0}')
%     title('Output SOP vs. \omega, for input SOP along \itx')
% end

% Compute spectrum of transmitted and received signals
Nblock = Npdg*Mct*(Nc+Npre);
deltaomegapdg = 2*pi*fsct/Nblock;
omegapdg = deltaomegapdg*fftshift([0:Nblock/2-1,-Nblock/2:-1]);

Setx = Nblock*fsct*fftshift(spectrum_periodogram(etx1,Nblock) + spectrum_periodogram(etx2,Nblock));
pdg_norm_fact = 2*pi*Ptx/deltaomegapdg/sum(Setx);
Setx = pdg_norm_fact*Setx;            % power spectral density of transmitted signal etx1 + etx2 (W/Hz)

Serx = Nblock*fsct*fftshift(spectrum_periodogram(erx1,Nblock) + spectrum_periodogram(erx2,Nblock));
Serx = pdg_norm_fact*Serx;            % power spectral density of received signal erx1 + erx2 (W/Hz)
                
% Plot spectrum of transmitted signal
gca = figure;
plot(1e-9/2/pi*omegapdg,10*log10(1e3*Setx));
h=get(gca,'children'); set(h,'linewidth',1.5)
title('Spectrum of transmitted signal');
xlabel('Frequency (GHz)');
ylabel('Power Spectral Density (dBm/Hz)');
gcax = get(gca,'CurrentAxes');
set(gcax,'XLim',[-fplotmax*1e-9,fplotmax*1e-9]);
text(omega(ucind(1))/2/pi*1e-9,max(10*log10(1e3*Setx)),'1')
text(omega(ucind(Nu/2))/2/pi*1e-9,max(10*log10(1e3*Setx)),num2str(Nu/2))
text(omega(ucind(Nu/2+1))/2/pi*1e-9,max(10*log10(1e3*Setx)),num2str(Nu/2+1))
text(omega(ucind(Nu))/2/pi*1e-9,max(10*log10(1e3*Setx)),num2str(Nu))

% Plot spectrum of received signal
gca = figure;
plot(1e-9/2/pi*omegapdg,10*log10(1e3*Serx));
h=get(gca,'children'); set(h,'linewidth',1.5)
title(['Spectrum of received signal at ' num2str(OSNR_dB_display,3) ' dB OSNR (2 pol. in ' num2str(osnrnoisebw/1e9,3) ' GHz BW)']);
xlabel('Frequency (GHz)');
ylabel('Power Spectral Density (dBm/Hz)');
gcax = get(gca,'CurrentAxes');
set(gcax,'XLim',[-fplotmax*1e-9,fplotmax*1e-9]);
text(omega(ucind(1))/2/pi*1e-9,max(10*log10(1e3*Serx)),'1')
text(omega(ucind(Nu/2))/2/pi*1e-9,max(10*log10(1e3*Serx)),num2str(Nu/2))
text(omega(ucind(Nu/2+1))/2/pi*1e-9,max(10*log10(1e3*Serx)),num2str(Nu/2+1))
text(omega(ucind(Nu))/2/pi*1e-9,max(10*log10(1e3*Serx)),num2str(Nu))

% Plot quantized encoded samples and D/A output (modulator drive signal)
figure
subplot(221)
plot(tsamp,vquant1i,'.',t,vmod1i,'-');
h=get(gca,'children'); set(h,'linewidth',1.5)
xlabel('Time (s)')
ylabel('x pol., in-phase (V)')
title('D/A: quant. input samples and interp. outputs')
subplot(222)
plot(tsamp,vquant1q,'.',t,vmod1q,'-');
xlabel('Time (s)')
ylabel('x pol., quadrature (V)')
subplot(223)
plot(tsamp,vquant2i,'.',t,vmod2i,'-');
xlabel('Time (s)')
ylabel('y pol., in-phase (V)')
subplot(224)
plot(tsamp,vquant2q,'.',t,vmod2q,'-');
xlabel('Time (s)')
ylabel('y pol., quadrature (V)')

% Plot filtered received signals and its samples
figure
subplot(221)
plot(tsamp,real(vsamp1),'.',t,real(vagc1),'-');
h=get(gca,'children'); set(h,'linewidth',1.5)
xlabel('Time (s)')
ylabel('x pol., in-phase (V)')
title(['A/D: filt. inputs and unquant. samples at ' num2str(OSNR_dB_display,3) ' dB OSNR (2 pol. in ' num2str(osnrnoisebw/1e9,3) ' GHz BW)']);
subplot(222)
plot(tsamp,imag(vsamp1),'.',t,imag(vagc1),'-');
xlabel('Time (s)')
ylabel('x pol., quadrature (V)')
subplot(223)
plot(tsamp,real(vsamp2),'.',t,real(vagc2),'-');
xlabel('Time (s)')
ylabel('y pol., in-phase (V)')
subplot(224)
plot(tsamp,imag(vsamp2),'.',t,imag(vagc2),'-');
xlabel('Time (s)')
ylabel('y pol., quadrature (V)')

% Plot samples of encoded signal without and with quantization
figure;
subplot(221)
plot(vdist1i,vdist1q,'.');
h=get(gca,'children'); set(h,'linewidth',1.5)
xlabel('x pol., in-phase (V)');
ylabel('x pol., quadrature (V)');
axis image
title('D/A: unquant. input samples')
subplot(222)
plot(vquant1i,vquant1q,'.');
h=get(gca,'children'); set(h,'linewidth',1.5)
xlabel('x pol., in-phase (V)');
ylabel('x pol., quadrature (V)');
axis image
title('D/A: quant. input samples')
subplot(223)
plot(vdist2i,vdist2q,'.');
h=get(gca,'children'); set(h,'linewidth',1.5)
xlabel('y pol., in-phase (V)');
ylabel('y pol., quadrature (V)');
axis image
subplot(224)
plot(vquant2i,vquant2q,'.');
h=get(gca,'children'); set(h,'linewidth',1.5)
xlabel('y pol., in-phase (V)');
ylabel('y pol., quadrature (V)');
axis image

% Plot samples of received signal without and with quantization
figure;
subplot(221)
plot(real(vsamp1),imag(vsamp1),'.');
h=get(gca,'children'); set(h,'linewidth',1.5)
xlabel('x pol., in-phase (V)');
ylabel('x pol., quadrature (V)');
axis image
title(['A/D: unquant. samples at ' num2str(OSNR_dB_display,3) ' dB OSNR']);
subplot(222)
plot(real(vquant1),imag(vquant1),'.');
h=get(gca,'children'); set(h,'linewidth',1.5)
xlabel('x pol., in-phase (V)');
ylabel('x pol., quadrature (V)');
axis image
title(['A/D: quant. samples at ' num2str(OSNR_dB_display,3) ' dB OSNR']);
subplot(223)
plot(real(vsamp2),imag(vsamp2),'.');
h=get(gca,'children'); set(h,'linewidth',1.5)
xlabel('y pol., in-phase (V)');
ylabel('y pol., quadrature (V)');
axis image
subplot(224)
plot(real(vquant2),imag(vquant2),'.');
h=get(gca,'children'); set(h,'linewidth',1.5)
xlabel('y pol., in-phase (V)');
ylabel('y pol., quadrature (V)');
axis image

% Plot step size and squared error
figure
subplot(211)
semilogy(1:Nsymb,muvec,'-',startdd*[1,1],[0.5*min(muvec),2*max(muvec)],'--');
h=get(gca,'children'); set(h,'linewidth',1.5)
title('Equalizer step size parameter vs. time');
xlabel('Symbol interval \itk'); ylabel('\mu(\itk\rm)');
text(startdd+4,0.5*max(muvec),['Start decision-directed at \itk\rm = ' num2str(startdd,3)])
legend('Step size parameter \mu(\itk\rm)')
subplot(212)
epsilonsq = reshape(sum(abs(epsilon).^2),Nu,Nsymb);
epsilonsqsm = [];
for k = 1:length(ch_sel)
    epsilonsqsm = [epsilonsqsm, smooth(epsilonsq(ch_sel(k),:),smoothlength)];
end
semilogy(1:Nsymb,epsilonsqsm');
h=get(gca,'children'); set(h,'linewidth',1.5)
xlabel('Symbol interval \itk'); ylabel('\mu(\itk\rm)');
ylabel('|\epsilon|^2(\itk\rm)');
legend(num2str(ch_sel',3))
title(['Squared error vs. time for selected subcarriers at ' num2str(OSNR_dB_adapt,3) ' dB OSNR (2 pol. in ' num2str(osnrnoisebw/1e9,3) ' GHz BW)']);

% Plot equalizer coefficient magnitudes
figure;

subplot(221)
plot(omegaucsort/2/pi*1e-9,20*log10(abs(reshape(Heq(1,1,ucindsort),1,Nu))),'.',...
    omegaucsort/2/pi*1e-9,20*log10(abs(reshape(Heq_ad(1,1,ucindsort),1,Nu))),'o');
h=get(gca,'children'); set(h,'linewidth',1.5)
xlabel('Subcarrier Frequency (GHz)');
ylabel('10 log_{10}(|\itH\rm_{eq, 11}|^{2}) (dB)')
legend('Fixed','Adaptive')
title('Magnitudes of one-tap equalizer coefficients')
axis([-fplotmax*1e-9 fplotmax*1e-9 eqmagmindB eqmagmaxdB]);
text(omega(ucind(1))/2/pi*1e-9,0.9*eqmagmaxdB,'1')
text(omega(ucind(Nu/2))/2/pi*1e-9,0.9*eqmagmaxdB,num2str(Nu/2))
text(omega(ucind(Nu/2+1))/2/pi*1e-9,0.9*eqmagmaxdB,num2str(Nu/2+1))
text(omega(ucind(Nu))/2/pi*1e-9,0.9*eqmagmaxdB,num2str(Nu))
subplot(222)
plot(omegaucsort/2/pi*1e-9,20*log10(abs(reshape(Heq(1,2,ucindsort),1,Nu))),'.',...
    omegaucsort/2/pi*1e-9,20*log10(abs(reshape(Heq_ad(1,2,ucindsort),1,Nu))),'o');
h=get(gca,'children'); set(h,'linewidth',1.5)
h=get(gca,'children'); set(h,'linewidth',1.5)
xlabel('Subcarrier Frequency (GHz)');
ylabel('10 log_{10}(|\itH\rm_{eq, 12}|^{2}) (dB)')
legend('Fixed','Adaptive')
axis([-fplotmax*1e-9 fplotmax*1e-9 eqmagmindB eqmagmaxdB]);
text(omega(ucind(1))/2/pi*1e-9,0.9*eqmagmaxdB,'1')
text(omega(ucind(Nu/2))/2/pi*1e-9,0.9*eqmagmaxdB,num2str(Nu/2))
text(omega(ucind(Nu/2+1))/2/pi*1e-9,0.9*eqmagmaxdB,num2str(Nu/2+1))
text(omega(ucind(Nu))/2/pi*1e-9,0.9*eqmagmaxdB,num2str(Nu))
subplot(223)
plot(omegaucsort/2/pi*1e-9,20*log10(abs(reshape(Heq(2,1,ucindsort),1,Nu))),'.',...
    omegaucsort/2/pi*1e-9,20*log10(abs(reshape(Heq_ad(2,1,ucindsort),1,Nu))),'o');
h=get(gca,'children'); set(h,'linewidth',1.5)
xlabel('Subcarrier Frequency (GHz)');
ylabel('10 log_{10}(|\itH\rm_{eq, 21}|^{2}) (dB)')
legend('Fixed','Adaptive')
axis([-fplotmax*1e-9 fplotmax*1e-9 eqmagmindB eqmagmaxdB]);
text(omega(ucind(1))/2/pi*1e-9,0.9*eqmagmaxdB,'1')
text(omega(ucind(Nu/2))/2/pi*1e-9,0.9*eqmagmaxdB,num2str(Nu/2))
text(omega(ucind(Nu/2+1))/2/pi*1e-9,0.9*eqmagmaxdB,num2str(Nu/2+1))
text(omega(ucind(Nu))/2/pi*1e-9,0.9*eqmagmaxdB,num2str(Nu))
subplot(224)
plot(omegaucsort/2/pi*1e-9,20*log10(abs(reshape(Heq(2,2,ucindsort),1,Nu))),'.',...
    omegaucsort/2/pi*1e-9,20*log10(abs(reshape(Heq_ad(2,2,ucindsort),1,Nu))),'o');
h=get(gca,'children'); set(h,'linewidth',1.5)
xlabel('Subcarrier Frequency (GHz)');
ylabel('10 log_{10}(|\itH\rm_{eq, 22}|^{2}) (dB)')
legend('Fixed','Adaptive')
axis([-fplotmax*1e-9 fplotmax*1e-9 eqmagmindB eqmagmaxdB]);
text(omega(ucind(1))/2/pi*1e-9,0.9*eqmagmaxdB,'1')
text(omega(ucind(Nu/2))/2/pi*1e-9,0.9*eqmagmaxdB,num2str(Nu/2))
text(omega(ucind(Nu/2+1))/2/pi*1e-9,0.9*eqmagmaxdB,num2str(Nu/2+1))
text(omega(ucind(Nu))/2/pi*1e-9,0.9*eqmagmaxdB,num2str(Nu))

% Plot equalizer coefficient phases
figure;
subplot(221)
plot(omegaucsort/2/pi*1e-9,phase(reshape(Heq(1,1,ucindsort),1,Nu)),'.',...
    omegaucsort/2/pi*1e-9,phase(reshape(Heq_ad(1,1,ucindsort),1,Nu)),'o');
h=get(gca,'children'); set(h,'linewidth',1.5)
xlabel('Subcarrier Frequency (GHz)');
ylabel('\angle(\itH\rm_{eq, 11}) (rad)')
legend('Fixed','Adaptive')
title('Phases of one-tap equalizer coefficients')
axis([-fplotmax*1e-9 fplotmax*1e-9 eqphsmin eqphsmax])
text(omega(ucind(1))/2/pi*1e-9,0.9*eqphsmax,'1')
text(omega(ucind(Nu/2))/2/pi*1e-9,0.9*eqphsmax,num2str(Nu/2))
text(omega(ucind(Nu/2+1))/2/pi*1e-9,0.9*eqphsmax,num2str(Nu/2+1))
text(omega(ucind(Nu))/2/pi*1e-9,0.9*eqphsmax,num2str(Nu))
subplot(222)
plot(omegaucsort/2/pi*1e-9,phase(reshape(Heq(1,2,ucindsort),1,Nu)),'.',...
    omegaucsort/2/pi*1e-9,phase(reshape(Heq_ad(1,2,ucindsort),1,Nu)),'o');
h=get(gca,'children'); set(h,'linewidth',1.5)
h=get(gca,'children'); set(h,'linewidth',1.5)
xlabel('Subcarrier Frequency (GHz)');
ylabel('\angle(\itH\rm_{eq, 12}) (rad)')
legend('Fixed','Adaptive')
axis([-fplotmax*1e-9 fplotmax*1e-9 eqphsmin eqphsmax])
text(omega(ucind(1))/2/pi*1e-9,0.9*eqphsmax,'1')
text(omega(ucind(Nu/2))/2/pi*1e-9,0.9*eqphsmax,num2str(Nu/2))
text(omega(ucind(Nu/2+1))/2/pi*1e-9,0.9*eqphsmax,num2str(Nu/2+1))
text(omega(ucind(Nu))/2/pi*1e-9,0.9*eqphsmax,num2str(Nu))
subplot(223)
plot(omegaucsort/2/pi*1e-9,phase(reshape(Heq(2,1,ucindsort),1,Nu)),'.',...
    omegaucsort/2/pi*1e-9,phase(reshape(Heq_ad(2,1,ucindsort),1,Nu)),'o');
h=get(gca,'children'); set(h,'linewidth',1.5)
xlabel('Subcarrier Frequency (GHz)');
ylabel('\angle(\itH\rm_{eq, 21}) (rad)')
legend('Fixed','Adaptive')
axis([-fplotmax*1e-9 fplotmax*1e-9 eqphsmin eqphsmax])
text(omega(ucind(1))/2/pi*1e-9,0.9*eqphsmax,'1')
text(omega(ucind(Nu/2))/2/pi*1e-9,0.9*eqphsmax,num2str(Nu/2))
text(omega(ucind(Nu/2+1))/2/pi*1e-9,0.9*eqphsmax,num2str(Nu/2+1))
text(omega(ucind(Nu))/2/pi*1e-9,0.9*eqphsmax,num2str(Nu))
subplot(224)
plot(omegaucsort/2/pi*1e-9,phase(reshape(Heq(2,2,ucindsort),1,Nu)),'.',...
    omegaucsort/2/pi*1e-9,phase(reshape(Heq_ad(2,2,ucindsort),1,Nu)),'o');
h=get(gca,'children'); set(h,'linewidth',1.5)
xlabel('Subcarrier Frequency (GHz)');
ylabel('\angle(\itH\rm_{eq, 22}) (rad)')
legend('Fixed','Adaptive')
axis([-fplotmax*1e-9 fplotmax*1e-9 eqphsmin eqphsmax])
text(omega(ucind(1))/2/pi*1e-9,0.9*eqphsmax,'1')
text(omega(ucind(Nu/2))/2/pi*1e-9,0.9*eqphsmax,num2str(Nu/2))
text(omega(ucind(Nu/2+1))/2/pi*1e-9,0.9*eqphsmax,num2str(Nu/2+1))
text(omega(ucind(Nu))/2/pi*1e-9,0.9*eqphsmax,num2str(Nu))

% Plot constellation diagrams of selected subcarriers
for s=1:length(ch_sel)
    figure;
    subplot(221)
    plot(real(Z1(ch_sel(s),:)),imag(Z1(ch_sel(s),:)),'.');
    h=get(gca,'children'); set(h,'linewidth',1.5)
    grid
    xlabel('x pol., in-phase (V)');
    ylabel('x pol., quadrature (V)');
    legend('Fixed Equalizer')
    axis image
    title(['Subcarrier ' num2str(ch_sel(s),3) ' at one-tap equalizer output at ' num2str(OSNR_dB_scan(end),3) ' dB OSNR (2 pol. in ' num2str(osnrnoisebw/1e9,3) ' GHz BW)']);
    subplot(223)
    plot(real(Z2(ch_sel(s),:)),imag(Z2(ch_sel(s),:)),'.');
    h=get(gca,'children'); set(h,'linewidth',1.5)
    grid
    xlabel('y pol., in-phase (V)');
    ylabel('y pol., quadrature (V)');
    legend('Fixed Equalizer')
    axis image
    subplot(222)
    plot(real(Z1ad(ch_sel(s),:)),imag(Z1ad(ch_sel(s),:)),'.');
    h=get(gca,'children'); set(h,'linewidth',1.5)
    grid
    xlabel('x pol., in-phase (V)');
    ylabel('x pol., quadrature (V)');
    legend('Adaptive Equalizer')
    axis image
    subplot(224)
    plot(real(Z2ad(ch_sel(s),:)),imag(Z2ad(ch_sel(s),:)),'.');
    h=get(gca,'children'); set(h,'linewidth',1.5)
    grid
    xlabel('y pol., in-phase (V)');
    ylabel('y pol., quadrature (V)');
    legend('Adaptive Equalizer')
    axis image
end;

% Plot bit-error probability vs. transmitted power
gca = figure;
semilogy(PtxdBm_ber,ber_scan,'-',PtxdBm_ber,ber_scan_ad,'--');
h=get(gca,'children'); set(h,'linewidth',1.5)
grid
xlabel('Transmitted power (dBm)');
ylabel('Bit-error probability');
legend('Fixed Equalizer','Adaptive Equalizer')
gcax = get(gca,'CurrentAxes');
title('Bit-error probability vs. transmitted power');

% Plot bit-error probability vs. OSNR
gca = figure;
semilogy(OSNR_dB_ber,ber_scan,'-',OSNR_dB_ber,ber_scan_ad,'--');
h=get(gca,'children'); set(h,'linewidth',1.5)
grid
xlabel('Optical Signal-to-Noise Ratio (dB)');
ylabel('Bit-error probability');
legend('Fixed Equalizer','Adaptive Equalizer')
gcax = get(gca,'CurrentAxes');
titlestr = ['Bit-error probability vs. OSNR in 2 polarizations in ' num2str(osnrnoisebw/1e9,3) ' GHz BW'];
title(titlestr);

toc