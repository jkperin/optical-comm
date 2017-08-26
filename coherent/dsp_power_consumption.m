%% Calculate power consumption
clear, clc, close all

addpath f

Rb = 224e9; % bit rate
process = 28e-9; % CMOS process
nb = 6; % resolution used in DSP
V = 0.8; % CMOS voltage
M = 4; % QAM order
Rs = Rb/(2*log2(M));
nadc = 6;
ros = 5/4;
Fs = Rs*ros;
Fa = 2.5e-12;
Ncd = 512;
Nfft = 4*Ncd;
Npmd = 7;
Bpmd = 128;
Btr = 256;
Bcr = 512;
Nf = 32;

P = PowerConsumption(Rb, M, process, nb, V);

Pdsp = Rb*[P.Eadc(Fa, nadc, Fs),...
    + 0*P.Ecd(Ncd, ros, Nfft),...
    + 0.5*P.Epmd(Npmd, Bpmd),...
    + P.Etr(Btr),...
    + P.Ecr(Nf, Bcr)]
    

sum(Pdsp)
