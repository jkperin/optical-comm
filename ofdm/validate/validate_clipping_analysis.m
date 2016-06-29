clear, clc, close all

addpath f   % functions path

rng('default') % initiate default random number generator

% Anonymous functions
fr = @(r) (1 + r.^2.*qfunc(r)).*(1 - qfunc(r)) + r.*exp(-r.^2/2)/sqrt(2*pi).*(2*qfunc(r)-1) - exp(-r.^2)/(2*pi);
Gr = @(r) (fr(r) - (1 - qfunc(r)).^2);

CLIP_RATIO = 1:0.2:3;

% Simulation parameters
sim.Pb = 1e-6;                     % target BER
sim.Nsymb = 2^15;                   % number of OFDM symbols
sim.Mct = 2;                        % oversampling rate for emulation of continuous time
NEP = 50e-12;                       % Noise equivalent power of the TIA at the receiver (A/sqrt(Hz))

%% AWGN Channel (without bandwidth limitation)
% Modulation parameters
ofdm.Nc = 64;                    % total number of subcarriers = FFT size (must be even)
ofdm.Nu = 52;                    % number of nonzero subcarriers (including complex conjugate)
ofdm.CS = 16;                     % Desired spectral efficiency 
ofdm.Rb = 100e9;                 % bit rate (total over two polarizations) (b/s)

ofdm.Ms = ofdm.Nc/ofdm.Nu;              % oversampling ratio
ofdm.Rs = 2*ofdm.Rb/log2(ofdm.CS);      % sampling rate without oversampling (Hz)

% Transmitter parameters 
tx.kappa = 1;                                              % dc slope
tx.Hl = @(f) ones(size(f));

ofdm.Npre_os = 0; 

% Receiver parameters
rx.R = 1;                      % responsivity
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

% snr to achieve target BER; this is assuming that clipping is negligible
snr = fzero(@(x) berqam(ofdm.CS, x) - sim.Pb, 20); 
snrl = 10^(snr/10);

Pnrx = snrl*ofdm.fs*rx.St/ofdm.Nc;  % Calculate equivalent power at each subcarrier at the receiver assuming only noise

for k = 1:length(CLIP_RATIO)   
    ofdm.clip_ratio = CLIP_RATIO(k);
    r = ofdm.clip_ratio;
      
    ofdm.Pn = (1/(1-qfunc(r)))^2*(Pnrx./(abs(rx.R*tx.kappa.*tx.Hl(ofdm.fc)).^2));
    tx.P = tx.kappa*tx.Hl(0)*(r*(1-qfunc(r)) + 1/sqrt(2*pi)*exp(-r^2/2))*sqrt(2*sum(ofdm.Pn));   
    
    [ber, ~, VarUn] = dc_ofdm(ofdm, tx, rx, sim, false);
    VarUa(k) = mean(VarUn);

    CNRawgn (k) = 10*log10(Pnrx*Gr(r)/(rx.St*ofdm.fs/ofdm.Nc));
    
    berawgn.count(k) = ber.count;
    berawgn.est(k) = ber.est;
end

clear ofdm tx ber

%% Fixed bit loading with bandwidth limitation and predistortion at tx
% Modulation parameters
ofdm.Nc = 64;                    % total number of subcarriers = FFT size (must be even)
ofdm.Nu = 52;                    % number of nonzero subcarriers (including complex conjugate)
ofdm.CS = 16;                     % Desired spectral efficiency 
ofdm.Rb = 100e9;                 % bit rate (total over two polarizations) (b/s)

ofdm.Ms = ofdm.Nc/ofdm.Nu;       % oversampling ratio
ofdm.Rs = 2*ofdm.Rb/log2(ofdm.CS);     % sampling rate without oversampling (Hz)

% Transmitter parameters 
tx.kappa = 1;                                              % dc slope
tx.fnl = 0.75*ofdm.Ms*ofdm.Rs/2;                                % cutoff frequency of the laser
tx.Hl = @(f) abs(1./(1 + 2*1j*f/tx.fnl - (f/tx.fnl).^2));  % laser freq. resp. (unitless) f is frequency vector (Hz)
                                                            % 2nd-order filter with damping = 1

ofdm.Npre_os = cyclic_prefix(tx, ofdm);  
                                                            
% Receiver parameters
rx.R = 1;                      % responsivity

% Time and frequency scales
ofdm.fs = ofdm.Rs*(ofdm.Nc + ofdm.Npre_os)/ofdm.Nu;                % chip rate (Hz)
fsct = sim.Mct*ofdm.fs;                                            % sampling frequency to emulate continuous time (Hz)
dt = 1/fsct;                                                       % time increment in emulating continuous time (s)
Ntot = sim.Mct*sim.Nsymb*(ofdm.Nc + ofdm.Npre_os);                 % total number of points simulated in continuous time
df = fsct/Ntot;                                                    % frequency increment in continuous time (Hz)   

sim.t = dt*(0:Ntot-1);                                             % continuous time scale
sim.f = -fsct/2:df:fsct/2-df;                                      % frequency range in continuous time                    
ofdm.fc = ofdm.fs/ofdm.Nc*(1:ofdm.Nu/2);                           % frequency at which subcarriers are located 


% snr to achieve target BER; this is assuming that clipping is negligible
snr = fzero(@(x) berqam(ofdm.CS, x) - sim.Pb, 20); 
snrl = 10^(snr/10);

Pnrx = snrl*ofdm.fs*rx.St/ofdm.Nc;  % Calculate equivalent power at each subcarrier at the receiver assuming only noise

for k = 1:length(CLIP_RATIO)   
    ofdm.clip_ratio = CLIP_RATIO(k);
    r = ofdm.clip_ratio;
       
    ofdm.Pn = (1/(1-qfunc(r)))^2*(Pnrx./(abs(rx.R*tx.kappa.*tx.Hl(ofdm.fc)).^2));
    tx.P = tx.kappa*tx.Hl(0)*(r*(1-qfunc(r)) + 1/sqrt(2*pi)*exp(-r^2/2))*sqrt(2*sum(ofdm.Pn));   
    
    [ber, ~, VarUn] = dc_ofdm(ofdm, tx, rx, sim, false);
    VarUf(k) = mean(VarUn);

    CNRfbl(k) = 10*log10(Pnrx*Gr(r)/(rx.St*ofdm.fs/ofdm.Nc));
    
    berfbl.count(k) = ber.count;
    berfbl.est(k) = ber.est;
end

clear ofdm tx ber

%% Variable bit loading with bandwidth limitation
% Modulation parameters
ofdm.Nc = 64;                    % total number of subcarriers = FFT size (must be even)
ofdm.Nu = 52;                    % number of nonzero subcarriers (including complex conjugate)
ofdm.CS = 16;                     % Desired spectral efficiency 
ofdm.Rb = 100e9;                 % bit rate (total over two polarizations) (b/s)

ofdm.Ms = ofdm.Nc/ofdm.Nu;       % oversampling ratio
ofdm.Rs = 2*ofdm.Rb/log2(ofdm.CS);     % sampling rate without oversampling (Hz)

% Transmitter parameters 
tx.kappa = 1;                                              % dc slope
tx.fnl = 0.75*ofdm.Ms*ofdm.Rs/2;                                % cutoff frequency of the laser
tx.Hl = @(f) abs(1./(1 + 2*1j*f/tx.fnl - (f/tx.fnl).^2));  % laser freq. resp. (unitless) f is frequency vector (Hz)
                                                           % 2nd-order filter with damping = 1

% Compute cyclic prefix length (assumes 2nd order filter with unit damping)
ofdm.Npre_os = cyclic_prefix(tx, ofdm);

% Receiver parameters
rx.R = 1;                      % responsivity

% Time and frequency scales
ofdm.fs = ofdm.Rs*(ofdm.Nc + ofdm.Npre_os)/ofdm.Nu;                % chip rate (Hz)
fsct = sim.Mct*ofdm.fs;                                            % sampling frequency to emulate continuous time (Hz)
dt = 1/fsct;                                                       % time increment in emulating continuous time (s)
Ntot = sim.Mct*sim.Nsymb*(ofdm.Nc + ofdm.Npre_os);                 % total number of points simulated in continuous time
df = fsct/Ntot;                                                    % frequency increment in continuous time (Hz)   

sim.t = dt*(0:Ntot-1);                                             % continuous time scale
sim.f = -fsct/2:df:fsct/2-df;                                      % frequency range in continuous time                    
ofdm.fc = ofdm.fs/ofdm.Nc*(1:ofdm.Nu/2);                           % frequency at which subcarriers are located 


% Levin-Campello's algorithm parameters
beta = 1; % Information Granularity: is the smallest incremental unit of information that can be transmitted
B = round(ofdm.Rb*((ofdm.Nc + ofdm.Npre_os)*1/ofdm.fs));

for k = 1:length(CLIP_RATIO)   
    ofdm.clip_ratio = CLIP_RATIO(k);
    r = ofdm.clip_ratio;
            
    ofdm.Pn = zeros(1, ofdm.Nu/2);
    G = abs(tx.kappa*tx.Hl(ofdm.fc)*rx.R).^2;
    for kk = 1
        GNR = G.*(1-qfunc(r))^2./(G.*Gr(r).*sum(2*ofdm.Pn)/(ofdm.Nc-1) + ofdm.fs*rx.St/ofdm.Nc); % abs(Hl(fc)).^2/N0;
        [bn, ofdm.CS, ofdm.Pn] = Levin_Campello_MA(B, beta, sim.Pb, GNR);
    end

    tx.P = tx.kappa*tx.Hl(0)*(r*(1-qfunc(r)) + 1/sqrt(2*pi)*exp(-r^2/2))*sqrt(2*sum(ofdm.Pn));
    
    CNRvbl(k) = 10*log10(mean(Pnrx*Gr(r)/(rx.St*ofdm.fs/ofdm.Nc)));
    
    [ber, ~, VarUn] = dc_ofdm(ofdm, tx, rx, sim, false);
    VarUv(k) = mean(VarUn);
    
    bervbl.count(k) = ber.count;
    bervbl.est(k) = ber.est;
end

clear ofdm tx ber

%% Figures
figure
% plot(CLIP_RATIO, log10(berawgn.est), '-k', CLIP_RATIO, log10(berawgn.count), ':*k')
% hold on
plot(CLIP_RATIO, log10(berfbl.est), '-r', CLIP_RATIO, log10(berfbl.count), ':sr')
hold on
plot(CLIP_RATIO, log10(bervbl.est), '-b', CLIP_RATIO, log10(bervbl.count), ':sb')
xlabel('clipping ratio', 'FontSize', 12)
ylabel('log_{10}(BER)', 'FontSize', 12)
legend('Estimated (Constant BL)', 'Counted (Constant BL)', 'Estimated (Variable BL)', 'Counted (Variable BL)')

figure
plot(CLIP_RATIO, CNRawgn, 'k')
hold on
plot(CLIP_RATIO, CNRfbl, 'r')
plot(CLIP_RATIO, CNRvbl, 'b')
xlabel('clipping ratio')
ylabel('Clipping-to-noise ratio (dB)')
legend('AWGN', 'Fix BL', 'Var BL')
% set(gca, 'ytick', -15:0)

r = CLIP_RATIO;

figure
plot(CLIP_RATIO, VarUa, '*:g')
hold on
plot(CLIP_RATIO, VarUf, '*:r')
plot(CLIP_RATIO, VarUv, '*:b')
plot(CLIP_RATIO, Gr(r), '-k')
xlabel('clipping ratio')
ylabel('Average normalized clipping noise variance')
legend('AWGN', 'Fix BL', 'Var BL', 'Theory')


figure
plot(CLIP_RATIO, abs(VarUa-Gr(r))./Gr(r), '*:g')
hold on
plot(CLIP_RATIO, abs(VarUf-Gr(r))./Gr(r), '*:r')
plot(CLIP_RATIO, abs(VarUv-Gr(r))./Gr(r), '*:b')
xlabel('clipping ratio')
ylabel('Normalized E(Var(U_n)) error |E(Var(U_n))-G(r)|/G(r)')
legend('AWGN', 'Fix BL', 'Var BL')

