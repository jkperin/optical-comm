%% Estimate receiver sensitivity of M-PAM system in amplified and unamplified (PIN or APD) IM-DD links 
clear, clc

addpath ../f/
addpath ../apd/

sim.BERtarget = 1.8e-4;
sim.rexdB = -15;
Lkm = 0:10;
sim.N = 2^14;
sim.Rb = 56e9; 
sim.Mct = 10;
wavelength = 1380e-9; 
alpha = 0; % chirp parameter
RIN = -150; % dB/Hz
modBW = (10:5:50)*1e9; % 30e9;
N0 = (20e-12)^2;
eq.type ='fixed td-sr-le';
eq.Ntaps = 9;

%% M-PAM
% PAM(M, bit rate, leve spacing : {'equally-spaced', 'optimized'}, pulse
% shape: struct containing properties of pulse shape 
mpam = PAM(2, sim.Rb, 'equally-spaced', select_pulse_shape('rect', 1));

sim.fs = mpam.Rs*sim.Mct;
[f, t] = freq_time(sim.N, sim.fs); 
sim.f = f;

Fiber = fiber(0);

% DAC
% Nhold = sim.Mct;
% hZOH = 1/Nhold*ones(1, Nhold);
filt = design_filter('bessel', 5, 0.7*mpam.Rs/(sim.fs/2));
% Hdac = filt.H(sim.f/sim.fs).*freqz(hZOH, 1, sim.f, sim.fs).*exp(1j*2*pi*sim.f/sim.fs*(Nhold-1)/2);
Hdac = filt.H(sim.f/sim.fs);

% Modulator
modfc = modBW/sqrt(sqrt(2)-1); % converts to relaxation frequency
Hrxfilt = design_filter('butter', 5, 0.7*mpam.Rs/(sim.fs/2));
Haa = Hrxfilt.H(sim.f/sim.fs);

varRIN = @(Plevel, noiseBW) 10^(RIN/10)*Plevel.^2*noiseBW;

% Photodiode
PD = pin(1, 10e-9);

% Optical amplifier
OptAmp = OpticalAmplifier('ConstantGain', 20, 5, wavelength); % Can't be ConstantOutputPower here
BWopt = sim.fs;
Npol =2;

PrecUnamp = zeros(size(modBW));
PrecAPDopt = zeros(size(modBW));
PrecAmp = zeros(size(modBW));
OSNRdB = zeros(size(modBW));
for k = 1:length(modBW)
    fprintf('modBW = %d GHz\n', modBW(k)/1e9);
    
    Hmod = 1./(1 + 2*1j*f/modfc(k) - (f/modfc(k)).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
    
    Fiber.L = 0;
    
    Hch = Hdac.*Hmod.*Fiber.Himdd(f, wavelength, alpha, 'small signal').*Haa;
    
    [~, eq] = equalize(eq, [], Hch, mpam, sim, false);
    
    Hrx = Haa.*eq.Hrx.*eq.Hff(sim.f/(mpam.Rs)); % antialiasing + matched + linear equalizer filters
    noiseBW = 1/2*trapz(f, abs(Hrx).^2); % filter noise BW (includes noise enhancement penalty)
    
    noiseSTDUnamp = PD.stdNoise(Haa.*eq.Hrx, eq.Hff(sim.f/(mpam.Rs)), N0, RIN, sim); % thermal + shot + RIN
    
    noiseSTDAmp = @(Plevel) sqrt(noiseBW*N0 + PD.varShot(Plevel, noiseBW) + varRIN(Plevel, noiseBW)... % thermal + shot + RIN
                    + PD.R^2*(OptAmp.varNoiseDD(Plevel/(OptAmp.Gain), noiseBW, BWopt, Npol))); % sig-spont + spont-spont

    PrecUnamp(k) = mpam.sensitivity(noiseSTDUnamp, sim.rexdB, sim.BERtarget);
    PrecAmp(k) = mpam.sensitivity(noiseSTDAmp, sim.rexdB, sim.BERtarget);
    OSNRdB(k) = 10*log10(dBm2Watt(PrecAmp(k))/(2*OptAmp.Ssp*OptAmp.BWref));                   
end

PrecAmp = PrecAmp - OptAmp.GaindB; % Refer required power to before the optical amplifier
% D = 1e6*abs(Fiber.D(wavelength)*Lkm);

figure(1), hold on, box on
hplot(1) = plot(modBW/1e9, PrecUnamp, '-', 'LineWidth', 2, 'DisplayName', sprintf('%d-PAM unamplified', mpam.M));
hplot(2) = plot(modBW/1e9, PrecAmp, '--', 'Color', get(hplot(1), 'Color'), 'LineWidth', 2, 'DisplayName', sprintf('%d-PAM amplified', mpam.M));
legend('-dynamiclegend')
xlabel('Modulator bandwidth (GHz)', 'FontSize', 12)
ylabel('Receiver sensitivity (dBm)', 'FontSize', 12)
set(gca, 'FontSize', 12)

figure(2), hold on, box on
plot(modBW/1e9, OSNRdB, 'Color', get(hplot(2), 'Color'), 'LineWidth', 2, 'DisplayName', sprintf('%d-PAM', mpam.M))
legend('-dynamiclegend')
xlabel('Modulator bandwidth (GHz)', 'FontSize', 12)
ylabel('OSNR required (dB)', 'FontSize', 12)
set(gca, 'FontSize', 12)

M = 4;
BERpam = @(P) 2*(M-1)/M*qfunc(P*sqrt(log2(M)/(M-1)^2*1/(sim.Rb*N0)));
PrxdBm = fzero(@(P) log10(BERpam(dBm2Watt(P))) - log10(sim.BERtarget), -20)