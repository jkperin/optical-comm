%% CD and PMD power consumption
% Pillai, B. S. G., Sedighi, B., Guan, K., Anthapadmanabhan, N. P., Shieh, 
% W., Hinton, K. J., & Tucker, R. S. (2014). End-to-end energy modeling and 
% analysis of long-haul coherent transmission systems. Journal of Lightwave 
% Technology, 32(18), 3093–3111.
clear, clc, close all

addpath ../../f/

% 100G DP-QPSK for 2400 km
Rb = 100e9; % 2*112e9; % bit rate
Rs = 120e9/4; % symbol rate
R = 0.83; % code rate
ros = 2; % oversampling ratio
M = 4;
Npol = 2;
nadc = 6; % nominal ADC resolution
nb = nadc;
Lkm = 2400; % fiber length in km
Fiber = fiber(Lkm*1e3);
Fiber.PMD = true;
Ncd = 0.032*(Rs/1e9)^2*Lkm*17e-6*1e3 % Fiber.Ntaps(Rs, ros, 1550e-9)
meanDGD = 50e-12;
Npmd = 10; %ceil(2*meanDGD*Rs)   
pt = 40; % CMOS node (nm)
V = 0.8; % supply voltage

%% Equivalent energy per operation. Values in J
EopG = 0.69e-15*pt*V^2; % energy per gate operation 
EopR = 3.43e-15*pt*V^2; % energy per register read and write 
EopRO = 1.71e-15*nb*pt*V^2; % average energy per nb -input ROMread per bit
EopA = 2.57e-15*nb*pt*V^2; % average energy per nb-bit adder op 
EopM = 2.57e-15*nb^2*pt*V^2; % average energy per nb-bit multiplier op

Eop = [EopG, EopR, EopRO, EopA, EopM];

%% CD
Nfft = 4*1024; % number of samples per FFT/IFFT operation
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

Ecd = 1e12*ros*sum(NopCD.*Eop.*EopAdj)/(R*Nn*log2(M)) % energy per information bit


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

Epmd = 1e12*sum(Eop.*NopPMD)/(R*Bpmd*log2(M))