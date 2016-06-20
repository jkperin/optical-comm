%% Power consumption analysis
% Pillai, B. S. G., Sedighi, B., Guan, K., Anthapadmanabhan, N. P., Shieh, 
% W., Hinton, K. J., & Tucker, R. S. (2014). End-to-end energy modeling and 
% analysis of long-haul coherent transmission systems. Journal of Lightwave 
% Technology, 32(18), 3093–3111.
clear, clc, close all

addpath ../../f/

%% 100G DP-QPSK for 2400 km
Rb = 112e9; % 2*112e9; % bit rate
Rs = 120e9/4; % symbol rate
ros = 2; % oversampling ratio
M = 4;
Npol = 2;
nadc = 6; % nominal ADC resolution
nb = nadc;
Lkm = 2400; % fiber length in km
Fiber = fiber(Lkm*1e3);
Fiber.PMD = true;
Ncd = Fiber.Ntaps(Rs, ros, 1550e-9)
Npmd = 10
pt = 40; % CMOS node (nm)
V = 0.8; % supply voltage

%% 40Gbit/s DP-QPSK
% Rb = 40e9; % 2*112e9; % bit rate
% Rs = 11.5; % symbol rate
% ros = 2; % oversampling ratio
% M = 4;
% Npol = 2;
% nadc = 6; % nominal ADC resolution
% nb = nadc;
% Ncd = 152
% Npmd = 5
% pt = 90; % CMOS node (nm)
% V = 1; % supply voltage

Prx = dBm2Watt(15); % received power (approximately equal to LO power)

%% Equivalent energy per operation. Values in J
EopG = 0.69e-15*pt*V^2; % energy per gate operation 
EopR = 3.43e-15*pt*V^2; % energy per register read and write 
EopRO = 1.71e-15*nb*pt*V^2; % average energy per nb -input ROMread per bit
EopA = 2.57e-15*nb*pt*V^2; % average energy per nb-bit adder op 
EopM = 2.57e-15*nb^2*pt*V^2; % average energy per nb-bit multiplier op

Eop = [EopG, EopR, EopRO, EopA, EopM];

%% DAC
% FD = 1; % DAC figure of merit
% nd = 6; % DAC resolution
% dacFs = 56e9; % DAC sampling frequency 
% 
% Edac = 4*FD*nd*Fs/Br; % energy per information bit

%% Local oscillator
Elo = 2.5/Rb; 

%% Photodiodes
R = 1; % responsivity
Vbias = 3.3; % bias voltage

Epd = 8*R*Vbias*Prx/Rb;

%% TIA and ACG
Etia = 1.88/(Rb*log2(M)); % energy per information bit

%% ADC
Fa = 2.5e-12; %  J/conv-step % ADC figure of merit
Fs = Rs*ros; % ADC sampling rate

Eadc = 4*Fa*nadc*Fs/Rb; % energy per information bit

%% DSP
%% CD
Nfft = 512; % number of samples per FFT/IFFT operation
Nn = Nfft - Ncd + 1;

% Operations required assuming 3 multiplicatation per complex
% multiplication
NopCDROM = 2*Nfft*log2(Nfft)-3*Nfft+8;
NopCDAdd = 6*Nfft*log2(Nfft)-3*Nfft+8;
NopCDreg = (4*Nn + 2*Ncd - 2)*(nadc+2);
NopCDGate = (nadc+2)*Nn*(6*log2(Nn/128)-4);

% [EopG, EopR, EopRO, EopA, EopM];
NopCD = Npol*[NopCDGate, NopCDreg, NopCDROM, NopCDAdd, 0];

EopAdj = [1 1 (nadc+2)/nb (nadc+2)/nb ((nadc+2)/nb).^2]; % adjusts power consumption per operation since CD requires higher resolution

Ecd = ros*sum(NopCD.*Eop.*EopAdj)/(Nn*log2(M)); % energy per information bit

%% Timing recovery
Btr = 256; % number of samples used in timing recovery

% Operations required assuming 3 multiplicatation per complex
% multiplication
NopTRMult = 7*Btr + 2;
NopTRAdd = 28*Btr + 1;
NopTRRom = 4.5*nb*Btr;
NopTRReg = nb*(16*Btr + 6);
NopTRGate = 5*nb*Btr;

% [EopG, EopR, EopRO, EopA, EopM];
NopTR = Npol*[NopTRGate, NopTRReg, NopTRRom, NopTRAdd, NopTRMult];

Etr = ros*sum(Eop.*NopTR)/(Btr*log2(M));

%% PMD
% Npmd = 5; % number of taps in PMD filters
Bpmd = 256; % is the number of symbols input per clock cycle

% Operations required assuming 3 multiplicatation per complex
% multiplication
NopPMDMult = 6*Npmd*(Bpmd+1);
NopPMDAdd = 14*Npmd*(Bpmd + 1) - 2;
NopPMDRom = 4*nb*Bpmd;
NopPMDReg = 2*nb*(5*Npmd-1);
NopPMDGate = 0;

% % Operations required assuming 4 multiplicatation per complex
% % multiplication
% NopPMDMult = 8*Npmd*(Bpmd+1);
% NopPMDAdd = 8*Npmd*(Bpmd + 1) - 2;
% NopPMDRom = 4*nb*Bpmd;
% NopPMDReg = 2*nb*(5*Npmd-1);
% NopPMDGate = 0;

% [EopG, EopR, EopRO, EopA, EopM];
NopPMD = Npol*[NopPMDGate, NopPMDReg, NopPMDRom, NopPMDAdd, NopPMDMult];

Epmd = sum(Eop.*NopPMD)/(Bpmd*log2(M));

%% Carrier recovery
Bcr = 512; % is the number of symbols input per clock cycle
Nf = 32; % the number of symbols used for frequency estimation and coarse phase estimation for every block of BCR symbols.

% Operations required assuming 3 multiplicatation per complex
% multiplication
NopCRMult = 6*Bcr + 3*Nf;
NopCRAdd = 22*Bcr + 7*Nf - 2 - 4*Bcr/Nf;
NopCRRom = nb*(9*Bcr + 3 + 6*Bcr/Nf);
NopCRReg = nb*Bcr + 2*nb;
NopCRGate = 2*nb*Bcr;

% [EopG, EopR, EopRO, EopA, EopM];
NopCR = Npol*[NopCRGate, NopCRReg, NopCRRom, NopCRAdd, NopCRMult];

Ecr = sum(Eop.*NopCR)/(Bcr*log2(M));

%% FEC
Efec = 2*200e-3/Rb; % mW in 28 nm CMOS
% Note: 200mW is the target for 4-PAM systems. Factor of 2 appears because
% of two polarizations

%% Total
Edsp = Ecd + Etr + Epmd + Ecr;

[Ecd Etr Epmd Ecr]*Rb


etadc = 0.93; % transceiver module power conversion efficienc
Erxbit = 1/etadc*(Elo + Epd + Etia + Eadc + Edsp + Efec); % power consumption per information bit

EdspTotal = sum([Ecd Etr Epmd Ecr]*Rb)
Erx = Erxbit*Rb


