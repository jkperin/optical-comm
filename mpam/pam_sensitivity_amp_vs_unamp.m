%% Estimate BER
clear, clc

addpath ../f/
addpath ../apd/

sim.BERtarget = 1.8e-4;
sim.rexdB = -15;
Lkm = 0:30;
sim.N = 2^14;
sim.Rb = 25e9; % 112e9;
sim.Mct = 10;
wavelength = 1550e-9;  %1380e-9;
alpha = 1; % chirp parameter
RIN = -150; % dB/Hz
modBW = 15e9; % 30e9;
N0 = (30e-12)^2;
eq.type ='fixed td-sr-le';
eq.Ntaps = 9;

%% M-PAM
% PAM(M, bit rate, leve spacing : {'equally-spaced', 'optimized'}, pulse
% shape: struct containing properties of pulse shape 
mpam = PAM(4, sim.Rb, 'equally-spaced', select_pulse_shape('rect', 1));

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
Hmod = 1./(1 + 2*1j*f/modfc - (f/modfc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
Hrxfilt = design_filter('butter', 5, 0.7*mpam.Rs/(sim.fs/2));
Hrx = Hrxfilt.H(sim.f/sim.fs);

varRIN = @(Plevel, noiseBW) 10^(RIN/10)*Plevel.^2*noiseBW;

% Photodiode
PD = pin(1, 10e-9);

% Optical amplifier
OptAmp = OpticalAmplifier('ConstantGain', 20, 5, wavelength); % Can't be ConstantOutputPower here
BWopt = sim.fs;
Npol =2;

PrecUnamp = zeros(size(Lkm));
PrecAmp = zeros(size(Lkm));
OSNRdB = zeros(size(Lkm));
for k = 1:length(Lkm)
    fprintf('L = %d\n', Lkm(k));
    
    Fiber.L = 1e3*Lkm(k);
    
    Hch = Hdac.*Hmod.*Fiber.Himdd(f, wavelength, alpha, 'small signal').*Hrx;
    
    [~, eq] = equalize(eq, [], Hch, mpam, sim, ~false);
    
    Htot = Hrx.*eq.Hrx.*eq.Hff(sim.f/(mpam.Rs)); % antialiasing + matched + linear equalizer filters
    noiseBW = 1/2*trapz(f, abs(Htot).^2); % filter noise BW (includes noise enhancement penalty)
    
    noiseSTDUnamp = PD.stdNoise(Hrx.*eq.Hrx, eq.Hff(sim.f/(mpam.Rs)), N0, RIN, sim); % thermal + shot + RIN
    
    noiseSTDAmp = @(Plevel) sqrt(noiseBW*N0 + PD.varShot(Plevel, noiseBW) + varRIN(Plevel, noiseBW)... % thermal + shot + RIN
                    + PD.R^2*(OptAmp.varNoiseDD(Plevel/(OptAmp.Gain), noiseBW, BWopt, Npol))); % sig-spont + spont-spont

    PrecUnamp(k) = mpam.sensitivity(noiseSTDUnamp, sim.rexdB, sim.BERtarget);
    PrecAmp(k) = mpam.sensitivity(noiseSTDAmp, sim.rexdB, sim.BERtarget);
    OSNRdB(k) = 10*log10(dBm2Watt(PrecAmp(k))/(2*OptAmp.Ssp*OptAmp.BWref));
end

PrecAmp = PrecAmp - OptAmp.GaindB; % Refer required power to before the optical amplifier
D = 1e6*abs(Fiber.D(wavelength)*Lkm);

figure(1), hold on, box on
hplot(1) = plot(D, PrecUnamp, '-', 'LineWidth', 2, 'DisplayName', '4-PAM unamplified');
hplot(2) = plot(D, PrecAmp, '--', 'Color', get(hplot(1), 'Color'), 'LineWidth', 2, 'DisplayName', '4-PAM amplified');
legend('-dynamiclegend')
xlabel('Dispersion (ps/nm)', 'FontSize', 12)
ylabel('Receiver sensitivity (dBm)', 'FontSize', 12)
set(gca, 'FontSize', 12)

figure(2), hold on, box on
plot(D, OSNRdB, 'Color', get(hplot(2), 'Color'), 'LineWidth', 2, 'DisplayName', '4-PAM')
legend('-dynamiclegend')
xlabel('Dispersion (ps/nm)', 'FontSize', 12)
ylabel('OSNR required (dB)', 'FontSize', 12)
set(gca, 'FontSize', 12)

M = 4;
BERpam = @(P) 2*(M-1)/M*qfunc(P*sqrt(log2(M)/(M-1)^2*1/(sim.Rb*N0)));
PrxdBm = fzero(@(P) log10(BERpam(dBm2Watt(P))) - log10(sim.BERtarget), -20)