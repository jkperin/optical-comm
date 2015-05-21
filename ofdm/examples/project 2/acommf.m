% acommf.m   Joseph M. Kahn  7/23/13, updated 8/16/13.
% Estimates the maximum achievable bit rate using ACO-OFDM in VCSEL MMF link.
% References
%   Notes of 7/23/13.
%   D.J.F. Barros and J.M. Kahn, "Comparison of OFDM and OOK in DD MMF Links", JLT 2011.
%   J. Armstrong and A.J. Lowery, "Power efficient optical OFDM", Electron Lett. 2006.
%   J. Armstrong et al, "Performance of ACO-OFDM in AWGN for IM/DD System", Globecom 2006.

clear

% General constants
j = sqrt(-1);
k = 1.38e-23;               
q = 1.602e-19;             
h = 6.626e-34;
c = 2.998e8;
Tamb = 300;                                 % ambient temperature

% Frequency response calculations and plots
lw = 1.5;                                   % plot linewidth
deltaffreqresp = 100e6;                     % frequency spacing
fmaxfreqresp = 50e9;                        % maximum frequency 
ffreqresp = 0:deltaffreqresp:fmaxfreqresp;  % frequency scale 

% Laser (assumes two-pole response)
lambda = 850e-9;                            % wavelength
Pavgmax = 1e-3;                             % maximum allowed output power (W)
kappa = 0.1;                                % d.c. slope efficiency (W/A)
Rl = 25;                                    % laser resistance (Ohm) 
fnl = 15.4e9;                               % relaxation oscillation frequency (Hz)
zetal = 0.95;                               % damping constant (should not be exactly 1)
Hl = @(f) 1./(1+j*zetal*f/fnl-(f/fnl).^2);  % norm. current-to-intensity freq. resp. (unitless)
                                            % f is frequency vector (Hz)
fl = ffreqresp(sum((abs(Hl(ffreqresp)).^2)>=0.5)+1);
                                            % -3-dB cutoff frequency (Hz);    
                                     
% Fiber (assumes single-pole response)
lossdB = 3;                             % attenuation, coupling loss and link margin (optical dB)
ff = 25e9;                              % -3-dB cutoff frequency (Hz);

Hf = @(f) 10^(-lossdB/10)./(1+j*f/ff);  % norm. current-to-intensity freq. resp. (W/W)
                                        % f is frequency vector (Hz)
                                        
% Photodetector and receiver (assumes infinite bandwidth, white noise)
eta = 0.9;                              % quantum efficiency (W/W)
R = eta*lambda*q/h/c;                   % responsivity (A/W)
RF = 100;                               % feedback resistance (Ohm)
FndB = 6; Fn = 10^(FndB/10);            % electrical noise figure (electrical dB)
Sn = @(f) 2*k*Tamb/RF*Fn*ones(size(f)); % two-sided input-ref. thermal noise PSD (A^2/Hz)
                                        % f is frequency vector (Hz)
                                        
% Error correction
GammadB = 6; Gamma = 10^(GammadB/10);% Coding gap (electrical dB), should be > 0 dB.
                                        % 12.2 dB corresponds to uncoded QAM with Pe = 1e-12
pclip = 1e-5;                          % Clipping probability at VCSEL and at A/D converter.
                                        % Should be chosen consistently with error correction code.
                                        
% Power allocation
powalloc = 1;                           % conventional power allocation
%powalloc = 2;                           % power allocation taking account of VCSEL frequency response

% OFDM (assumes oversampling rate Ms = 1)
Rs = 64e9;                              % sampling rate (Hz)
N = 64;                                 % FFT length (should be a power of 2)
Deltaf = Rs/N;                          % subcarrier spacing (Hz)
m = 1:N/4;                              % indices of used positive-frequency subcarriers
fu = (2*m-1)*Rs/N;                      % frequencies of used positive-frequency subcarriers

% Analog-to-digital converter
b = 5;                          % number of quantizer bits
iclipu = 3.1e-3                         % user's guess of clipping level (A). Adjust until becomes consistent with program output.
Sq = @(f) 2^(-2*b)/12*iclipu^2/Rs*ones(size(f));
                                        % two-sided input-referred equivalent quantization noise PSD (A^2/Hz)
                                        % f is frequency vector (Hz)

% Compute cyclic prefix length (assumes laser has two poles and fiber has one pole)
% Does not include CP overhead or oversampling. 
frac_incl = 0.999;                      % fraction of energy to be included within CP length
tcp = 1/Rs*(0:1:100);                   % vector of sampling intervals
c1 = 2*pi*fnl*(-zetal+sqrt(zetal^2-1));
c2 = 2*pi*fnl*(-zetal-sqrt(zetal^2-1));
Ml = pi*fnl/sqrt(zetal^2-1);
hl = Ml*(exp(c1*tcp)-exp(c2*tcp));      % laser impulse response
hf = exp(-2*pi*ff*tcp);                 % fiber impulse response 
htot = conv(hl,hf);                     % total impulse response
en_frac = cumsum(htot.^2)/sum(htot.^2); % fraction of energy in impulse response vs. sampling interval
L = sum(en_frac<frac_incl)+1;           % number of samples that include the desired fraction of energy
                                        % minimum cyclic prefix length is L-1
CPpendB = 10*log10((N+L-1)/N);          % cyclic prefix penalty (electrical dB)
                                        
% Compute noise-to-gain ratio
NGR = (N+L-1)/N*Deltaf/(R*kappa/2)^2*...% noise-to-gain ratios in subcarriers
    (Sn(fu)+ Sq(fu))./abs(Hl(fu)).^2./abs(Hf(fu)).^2;

% Optimize power allocation

if powalloc == 1;                               % conventional power allocation
objfn = @(A) sqrt(kappa^2/pi*sum(abs(Hl(fu)).^2.*max(0,A-Gamma*NGR)))-Pavgmax;
                                                % objective function that should equal zero
Aopt = fzero(objfn,0);                          % optimized water level
Popt = max(0,Aopt-Gamma*NGR);                   % optimized electrical powers of subcarriers (A^2)

elseif powalloc == 2;                           % power allocation taking account of VCSEL freq. resp.
objfn = @(A) sqrt(kappa^2/pi*sum(abs(Hl(fu)).^2.*max(0,A./abs(Hl(fu)).^2-Gamma*NGR)))-Pavgmax;
                                                % objective function that should equal zero
Aopt = fzero(objfn,0);                          % optimized water level
Popt = max(0,Aopt./abs(Hl(fu)).^2-Gamma*NGR);   % optimized electrical powers of subcarriers (A^2)

end

% Compute clipping noise 
objfn = @(x) pclip/2-0.5*erfc(x/sqrt(2));
iclipp = fzero(objfn,0)*sqrt(2*kappa^2*R^2*sum(abs(Hl(fu)).^2.*abs(Hf(fu)).^2.*Popt))
                                                % clipping level (A). Adjust iclipu until it becomes consistent with this value.


% Compute SNRs
SNR = Popt./NGR;                        % electrical SNRs in subcarriers (linear units)

% Compute bit rates
Bn = log2(1+SNR/Gamma);                  % bits per symbol in subcarriers (bit)
Btot = Deltaf*sum(Bn);                   % total bit rate (bit/s)

% Compute VCSEL output power and power dissipation
Pavg = sqrt(kappa^2/pi*sum(abs(Hl(fu)).^2.*Popt));  % VCSEL optical output power (W)
objfn = @(x) pclip/2-0.5*erfc(x/sqrt(2));
Pclip = fzero(objfn,0)*sqrt(2*kappa^2*sum(abs(Hl(fu)).^2.*Popt));
                                                    % VCSEL clipping power, exceeded with prob. pclip (W) 
Pdiss = Rl*sum(Popt);                               % VCSEL electrical power dissipation (W)

% Plot results

% power allocation and system performance
figure(1)
clf

subplot(221)
plot(fu/1e9,1e3*sqrt(Popt),'b.')
xlabel('Subcarrier Frequency (GHz)')
ylabel('R.M.S. Current P_{n}^{1/2}(mA)');
title('VCSEL Drive Current Allocations');
axis([0 Rs/2/1e9 0 0.12*ceil(max(1e4*sqrt(Popt)))]);

subplot(222)
plot(fu/1e9, 10*log10(SNR),'b.')
xlabel('Subcarrier Frequency (GHz)')
ylabel('SNR_{n} (dB)');
title('Signal-to-Noise Ratios');
axis([0 Rs/2/1e9 0 0.12*ceil(max(100*log10(SNR)))]);

subplot(223)
plot(fu/1e9, b,'b.')
xlabel('Subcarrier Frequency (GHz)')
ylabel('Bits per Symbol B_{n}(bit)');
title('Bit Allocations');
axis([0 Rs/2/1e9 0 0.12*ceil(10*max(b))]);

subplot(224)
axis off
text(-0.1,1,['R_{s} = ' num2str(Rs/1e9,3) ' GHz, b = ' num2str(b,3) ' bits']);
text(-0.1,0.8,['N = ' num2str(N,3) ', L = ' num2str(L,3) ', CP penalty = ' num2str(CPpendB,2) ' dB']);
text(-0.1,0.6,['Power Allocation ' num2str(powalloc,1) ', \Gamma = ' num2str(10*log10(Gamma),3) ' dB, p_{clip} = ' num2str(pclip,3)]);
text(-0.1,0.4,['B_{tot} = ' num2str(Btot/1e9,3) ' Gbit/s']);
text(-0.1,0.2,['P_{avg} = ' num2str(Pavg/1e-3,3) ' mW, P_{clip} = ' num2str(Pclip/1e-3,3) ' mW, P_{diss} = ' num2str(Pdiss/1e-3,3) ' mW']);
text(-0.1,0.0,[ 'f_{VCSEL} = ' num2str(fl/1e9,3) ' GHz, f_{MMF} = ' num2str(ff/1e9,3) ' GHz, Loss = ' num2str(lossdB,3) ' dB']);

% laser and fiber
figure(2)
clf

subplot(221)
plot(ffreqresp/1e9,10*log10(abs(Hl(ffreqresp)).^2),'-b');
h=get(gca,'children'); set(h,'linewidth',lw)
xlabel('Frequency (GHz)')
ylabel('10 log_{10}(|H_{laser}(f)|^{2}) (dB)')
title('VCSEL Response')
axis([0 fmaxfreqresp/1e9 -20 10]);

subplot(222)
plot(ffreqresp/1e9,10*log10(abs(Hf(ffreqresp)).^2),'-b');
h=get(gca,'children'); set(h,'linewidth',lw)
xlabel('Frequency (GHz)')
ylabel('10 log_{10}(|H_{fiber}(f)|^{2}) (dB)')
title('Fiber Response (including loss)')
axis([0 fmaxfreqresp/1e9 -20 10]);                                        
