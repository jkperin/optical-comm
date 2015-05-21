%% BER floor
clear, clc, close all

addpath f   % functions path

rng('default') % initiate default random number generator

% Anonymous functions
fr = @(r) (1 + r.^2.*qfunc(r)).*(1 - qfunc(r)) + r.*exp(-r.^2/2)/sqrt(2*pi).*(2*qfunc(r)-1) - exp(-r.^2)/(2*pi);
Gr = @(r) (fr(r) - (1 - qfunc(r)).^2);

% --
% Simulation parameters
sim.Nsymb = 2^17;                   % number of OFDM symbols
sim.Mct = 1;                        % oversampling rate for emulation of continuous time
NEP = 0;                            % Noise equivalent power of the TIA at the receiver (A/sqrt(Hz))

% Modulation parameters
ofdm.ofdm = 'dc-ofdm';
ofdm.Nc = 64;                    % total number of subcarriers = FFT size (must be even)
ofdm.Nu = 52;                    % number of nonzero subcarriers (including complex conjugate)
ofdm.CS = 16;                    % Desired spectral efficiency 
ofdm.Rb = 100e9;                 % bit rate (total over two polarizations) (b/s)

ofdm.Ms = ofdm.Nc/ofdm.Nu;       % oversampling ratio
ofdm.Rs = 2*ofdm.Rb/log2(ofdm.CS);     % sampling rate without oversampling (Hz)

% Transmitter parameters 
tx.kappa = 1;                                       % dc slope
tx.fnl = 30e9;                                           % cutoff frequency of the laser
tx.Hl = @(f) abs(1./(1 + 2*1j*f/tx.fnl - (f/tx.fnl).^2));  % laser freq. resp. (unitless) f is frequency vector (Hz)
% tx.Hl = @(f) ones(size(f));

% Compute cyclic prefix length (assumes 2nd order filter with unit damping)
ofdm.Npre_os = cyclic_prefix(tx, ofdm);

% Receiver parameters
rx.R = 1;                       % responsivity

rx.St = rx.R^2*NEP^2/2;                      % two-sided psd of thermal noise at the receiver
    
% Time and frequency scales
ofdm.fs = ofdm.Rs*(ofdm.Nc + ofdm.Npre_os)/ofdm.Nu;                % chip rate (Hz)
fsct = sim.Mct*ofdm.fs;                                            % sampling frequency to emulate continuous time (Hz)
dt = 1/fsct;                                                       % time increment in emulating continuous time (s)
Ntot = sim.Mct*sim.Nsymb*(ofdm.Nc + ofdm.Npre_os);                 % total number of points simulated in continuous time
df = fsct/Ntot;                                                    % frequency increment in continuous time (Hz)   

sim.t = dt*(0:Ntot-1);                                             % continuous time scale
sim.f = -fsct/2:df:fsct/2-df;                                      % frequency range in continuous time                    
ofdm.fc = ofdm.fs/ofdm.Nc*(1:ofdm.Nu/2);                           % frequency at which subcarriers are located 
         
%% 
CLIP_RATIO = 1:0.2:2.2;
for k = 1:length(CLIP_RATIO)
    ofdm.clip_ratio = CLIP_RATIO(k);
    r = ofdm.clip_ratio;

    ofdm.Pn = 1e-3*ones(1, ofdm.Nu/2);

    tx.P = tx.kappa*tx.Hl(0)*(r*(1-qfunc(r)) + 1/sqrt(2*pi)*exp(-r^2/2))*sqrt(2*sum(ofdm.Pn));   

    [ber, ~, ~] = dc_ofdm(ofdm, tx, rx, sim, ~true);

    berest(k) = ber.est;
    bercount(k) = ber.count;
    interval(k,:) = ber.interval;
end
    

%% Figures
r = CLIP_RATIO;
snrfloor = 10*log10((1-qfunc(r)).^2./Gr(r));
berfloor16qam = berqam(16, snrfloor);
berfloor32qam = berqam(32, snrfloor);
berfloor64qam = berqam(64, snrfloor);

figure
plot(r, log10(berfloor16qam), 'k', r, log10(berfloor32qam), 'r', r, log10(berfloor64qam), 'b')
hold on
plot(r, log10(berest), '--*k', r, log10(bercount), '--dk')
plot(r, interval(:,1), 'xk', r, interval(:,2), 'xk')
axis([r(1) r(end) -6 0])
xlabel('clipping ratio', 'FontSize', 12)
ylabel('log_{10}(BER_{floor})', 'FontSize', 12)

legend('16-QAM', '32-QAM', '64-QAM')