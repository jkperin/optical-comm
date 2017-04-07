%% Calculate receiver sensitivity and OSNR required for SSB-OFDM
clear, clc
addpath f/
addpath ../f/
addpath ../apd/

sim.OFDM = 'SSB-OFDM'; % either 'DC-OFDM' or 'ACO-OFDM'
sim.Rb = 112e9; % Bit rate
sim.BERtarget = 1.8e-4; % target BER
sim.Mct = 5; % Oversampling ratio to emulate continuous time
Lkm = 0:15; % Fiber length in km

N0 = (30e-12)^2; % thermal noise PSD
wavelength = 1550e-9; 
modBW = 30e9; % modulator bandwidth

quantiz = true; % whether quantization is assumed

%% OFDM 
% OFDM constructor: ofdm(Nc, Nu, CS, Rb, power_allocation_type (optional, default = 'palloc')
% Nc : number of subcarriers
% Nu : number of used subcarriers
% CS : nominal constellation size
% Rb : bit rate (b/s)
% power_allocation_type : {'palloc', 'preemphasis'}
ofdm = ofdm(256, 208, 16, 112e9, 'SSB', 'Levin-Campello-MA'); 
ENOB = 5;
rclip = 3; % clipping ratio (controls the DC-bias of the DC-OFDM signal)

ofdm.set_cyclic_prefix(5, 5); % set cyclic prefix length. Should be consistent with channel memory length
noiseBW = ofdm.fs/2;

sim.fs = ofdm.fs*sim.Mct;

% Modulator
filt = design_filter('butter', 5, noiseBW/(sim.fs/2));
Nhold = sim.Mct;
hZOH = 1/Nhold*ones(1, Nhold);
Hdac = filt.H(ofdm.fc/sim.fs).*freqz(hZOH, 1, ofdm.fc, sim.fs).*exp(1j*2*pi*ofdm.fc/sim.fs*(Nhold-1)/2);

Mod.filt = design_filter('two-pole', modBW, sim.fs);
% Tx.Mod.filt = design_filter('butter', 5, Tx.Mod.BW/(sim.fs/2));
Hmod = Mod.filt.H(ofdm.fc/sim.fs);

% modfc = modBW/sqrt(sqrt(2)-1); % converts to relaxation frequency
% Hmod = 1./(1 + 2*1j*ofdm.fc/modfc - (ofdm.fc/modfc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)

% Fiber
Fiber = fiber(0);

% Optical amplifier
OptAmp = OpticalAmplifier('ConstantGain', 20, 5, wavelength);
% OptAmp = OpticalAmplifier('ConstantOutputPower', 0, 5, wavelength);
BWopt = sim.fs;
Npol =2;

% Photodiode
PD = pin(1, 10e-9);

% ADC
Hadc = filt.H(ofdm.fc/sim.fs);
     
%% Amplified system
% Note: Prx is referred to the input of the optical amplifier
% Channel frequency response
Hch = Hdac.*Hmod.*Hadc;
% Hch = ones(size(ofdm.fc));
gamma = 1e-7;
CSPRdB = 0:20;
Bref = 12.5e9;
for k = 1:length(CSPRdB)
    CSPR = 10^(CSPRdB(k)/10);
    tol = Inf;
    Ps = 1e-6;
    Pc = CSPR*Ps;
    n = 1;
    while tol > 1e-3 && n < 50        
        varNoiseAmp = 1/(ofdm.Nc)*(abs(Hadc).^2*4*PD.R^2*OptAmp.Gain*(Ps+Pc)*OptAmp.Ssp*ofdm.fs/2 ... % amp
            + 1/3*rclip^2*(2*PD.R^2*OptAmp.Gain^2*Pc*Ps)*2^(-2*ENOB)... % quantization
            + gamma*PD.R^2*OptAmp.Gain^2*Ps);
        G = PD.R*OptAmp.Gain*sqrt(Pc);
        ofdm.power_allocation(Hch*G, varNoiseAmp, sim.BERtarget, false);
        Pnrx = ofdm.Pn; % referred to the input of the amplifier
        Psnew = sum(Pnrx);
        tol = (Psnew-Ps)/Psnew;
        Ps = Psnew;
        Pc = Ps*CSPR;
        n = n + 1;
    end
    n
        
    % OSNR calculation: assumes that total power is equal to carrier power
    OSNRdBreq(k) = 10*log10(OptAmp.Gain*Pc/(2*OptAmp.Ssp*Bref));
end

figure(1), hold on, box on
plot(CSPRdB, OSNRdBreq)
xlabel('CSPR (dB)')
ylabel('OSNR (dB)')
