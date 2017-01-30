%% Theoretic BER curves
clear, clc, close all

addpath ../f/
addpath ../mpam/
addpath ../ofdm/
addpath ../ofdm/f/
addpath f/

%% Simulation parameters
sim.Rb = 56e9;
sim.Mct = 1;      
sim.BERtarget = 1.8e-4;
PrxdBm = 13; % expected received power in dBm
Prx = dBm2Watt(PrxdBm);

% Fiber
Fiber = fiber(0e3);

%% Amplifier
%% ======================== Optical Amplifier =============================
% Constructor: OpticalAmplifier(GaindB, NF, lamb, maxGain (optional))
% GaindB : Gain in dB
% NF : Noise figure in dB
% lamb : wavelength in m
% maxGain = maximum amplifier gain in dB, default is Inf
EDFA = OpticalAmplifier(20, 5, 1550e-9);
OptAmpOutPowerdBm = 0; % output power after amplifier
% Note: the amplifier here operates in the constant output power mode,
% where the output power after amplification is set to Rx.AmpOutPowerdBm

%% OFDM 
% OFDM constructor: ofdm(Nc, Nu, CS, Rb, power_allocation)
% Nc : number of subcarriers
% Nu : number of used subcarriers
% CS : nominal constellation size
% Rb : bit rate (b/s)
% power_allocation (optional, default = 'palloc') : type of power
% allocation {'preemphasis', 'palloc'}
ofdm = ofdm(256, 208, 16, sim.Rb, 'palloc'); 
Ncp = 12;
ofdm.set_cyclic_prefix(Ncp/2, Ncp/2);

%% Power allocation 
Hpreemph = 10.^(polyval([-0.0013 0.5846 1.5859], ofdm.fc/1e9)/20); 
Hmod = 1./Hpreemph; % MZM frequency response
Hmod = 1;

Himdd = Fiber.Himdd(ofdm.fc, 1550e-9, 0, 'large signal');

Hch = Hmod.*Himdd;
OSNRdB = 25;
BWref = 12.5e9; % reference bandwidth used to measure OSNR
OSNR = 10^(OSNRdB/10);
Ssp = Prx/(2*BWref*OSNR); % Amplifier one-sided ASE PSD per polarization 
varSigSpont = 4*Prx*Ssp*ofdm.fs/2;

ofdm.power_allocation(Hch, varSigSpont*ones(size(ofdm.fc))/ofdm.Nc, sim.BERtarget, true);

%% MPAM
mpam = PAM(4, sim.Rb);

OSNRdB = 20:40;

berPAM = pam_ber_from_osnr(mpam.M, OSNRdB, mpam.Rs/2);

ofdm.Pn = ofdm.Pn.*abs(Hch).^2; % channel effect
Pmean = ofdm.dc_bias(ofdm.Pn, 3);
ofdm.Pn = ofdm.Pn*(Prx/Pmean)^2;

[Prx, ofdm.dc_bias(ofdm.Pn, 3)] 

for k = 1:length(OSNRdB)
    OSNR = 10^(OSNRdB(k)/10);
    
    Ssp = Prx/(2*BWref*OSNR); % Amplifier one-sided ASE PSD per polarization 
    
    % Calculate unbiased estimated SNR
    SNRn = 10*log10(ofdm.Pn/(4*Prx*Ssp*ofdm.fs/(2*ofdm.Nc))); 

    berest = zeros(size(SNRn));
    for kk = 1:length(SNRn)
        if ofdm.CSn(kk) == 0
            continue
        else
            berest(kk) = berqam(ofdm.CSn(kk), SNRn(kk));
        end
    end

    berOFDM(k) = sum(berest.*ofdm.bn)/sum(ofdm.bn);
%     M = ofdm.CS;
%     bertheory(k) = 4/log2(M)*(1 - 1/sqrt(M))*qfunc(Prx/3*sqrt(3/(M-1)*log2(M)/(4*Prx*Ssp*sim.Rb)));
end


figure
semilogy(OSNRdB, berPAM, 'LineWidth', 4)
hold on
semilogy(OSNRdB, berOFDM, 'LineWidth', 4)
% plot(OSNRdB, log10(bertheory), 'LineWidth', 2)
xlabel('OSNR (dB)', 'FontSize', 18)
ylabel('BER', 'FontSize', 18)
legend(sprintf('%d-PAM', mpam.M), sprintf('%d-QAM OFDM', ofdm.CS))
% title(sprintf('Theoretical BER curves for %d-Gbit/s PAM and OFDM', sim.Rb/1e9))
set(gca, 'FontSize', 18)
set(gca, 'ytick', [1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1])
grid on
axis([OSNRdB([1 end]) 1e-8 1e-1])