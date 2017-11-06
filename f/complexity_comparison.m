clear, clc, close all


addpath ../mpam/
addpath ../ofdm/
addpath ../coherent/f/

wavelength = 1550e-9;
alpha = 0;
Fiber = fiber(4.7e3);
Rb = 112e9;

% L = 4.7e3/2;
L = 4.7e3;
% L = 80e3;
% L = (1:25)*1e3;

P6 = PowerConsumption(112e9, 4, 28e-9, 6)
Eavg6 = (P6.Eadd + P6.Emult)/2

P7 = PowerConsumption(112e9, 4, 28e-9, 7)
Eavg7 = (P7.Eadd + P7.Emult)/2

for l = 1:length(L)
    Fiber.L = L(l);
    D(l) = Fiber.L/1e3*Fiber.D(wavelength)*1e6;
    [OPbit, Fs] = operations_per_bit(Rb, Fiber, wavelength);
    
    N4PAMbitFD(l) = OPbit.N4PAMbitFD;
    N4PAMbitTD(l) = OPbit.N4PAMbitTD;
    NOFDMbit(l) = OPbit.NOFDMbit;
    NSSBOFDMbit(l) = OPbit.NSSBOFDMbit;
    N2DStokesbit(l) = OPbit.N2DStokesbit;
    N3DStokesbit(l) = OPbit.N3DStokesbit;
    NQPSKbit(l) = OPbit.NQPSKbit;
    N16QAMbit(l) = OPbit.N16QAMbit;
    NKK4PAMbitTD(l) = OPbit.NKK4PAMbitTD;
    NKK4PAMbitFD(l) = OPbit.NKK4PAMbitFD;
    
    
    figure, box on
    y = [min(N4PAMbitTD(l), NKK4PAMbitFD(l)) NOFDMbit(l) NSSBOFDMbit(l) min(NKK4PAMbitTD(l), OPbit.NKK4PAMbitFD)...
        N2DStokesbit(l), N3DStokesbit(l),...
        NQPSKbit(l), N16QAMbit(l)];
    barh(1:length(y), y)
    c = {'4-PAM','DC-OFDM','SSB-OFDM', 'KK 4-PAM',...
        'Stokes 2D', 'Stokes 3D',...
        'DP-QPSK', 'DP-16QAM'};
    set(gca,'yticklabel', c);
    xlabel('Number of real operations per bit')
    annotation(gcf,'textbox',...
    [0.626190476190476 0.156142857142857 0.256547619047619 0.106349206349206],...
    'LineStyle','none', 'String',{'Bit rate = 112 Gbit/s','D = 80 ps/nm'});

    % Order
    % 4-PAM, OFDM, SSB-OFDM, KK-PAM,...
    % Stokes 2D, Stokes 3D,
    % DP-QPSK, DP-16QAM    
    
    Eavg = [Eavg6, Eavg7, Eavg7, Eavg7,...
        Eavg7, Eavg7,...
        Eavg6, Eavg7]; 
    
    Resolution = [6 7 7 7,...
        7 7,...
        6 7];
    
    ResolutionDAC = [2 7 7 2,...
    2 2,...
    1 2]; % Only OFDM is assumed to need high-resolution DACs
    
    fs = [Fs.PAM4, Fs.DCOFDM, Fs.SSBOFDM, Fs.KKPAM4,...
        Fs.Stokes2D, Fs.Stokes3D,...
        Fs.QPSK, Fs.QAM16];
    
    NconvDAC = 1/4*[1 1 1 1,...
        2 3,...
        4 4];
    
    NconvADC = 1/4*[1 1 1 1,...
        4 4,...
        4 4];
            
    Pdac = P6.Rb*P6.Edac(1.56e-12, ResolutionDAC, fs).*NconvDAC;
    Padc = P6.Rb*P6.Eadc(2.5e-12, Resolution, fs).*NconvADC;
        
    Pconsumption = [Pdac; %DAC
        Padc;... % ADC
        Rb*Eavg.*y];
    
    figure, box on
    barh(1:length(y), Pconsumption.', 'stacked')
    set(gca, 'yticklabel', c);
    xlabel('Power consumption in 28 nm CMOS (W)')
    annotation(gcf,'textbox',...
    [0.626190476190476 0.156142857142857 0.256547619047619 0.106349206349206],...
    'LineStyle','none', 'String',{'Bit rate = 112 Gbit/s','D = 80 ps/nm'});
    legend('DAC', 'ADC', 'DSP')
    xlim([0 40])
    
end

Pconsumption