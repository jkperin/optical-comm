%% Estimate receiver sensitivity of M-PAM system in amplified and unamplified (PIN or APD) IM-DD links 
clear, clc

addpath ../f/
addpath ../apd/

sim.BERtarget =3.8e-3;
sim.rexdB = -15;
Lkm = 0:10;
sim.N = 2^14;
sim.Rb = 112e9; 
sim.Mct = 10;
wavelength = 1380e-9; 
alpha = 0; % chirp parameter
RIN = -150; % dB/Hz
modBW = 30e9; % 30e9;
N0 = (20e-12)^2;
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
% constructor: apd(GaindB, ka, BW, R, Id)
APD = apd(10, 0.18, [24e9 290e9], 0.74, 40e-9);
% APD = apd(10, 0.18, Inf, 0.74, 40e-9);

% Optical amplifier
OptAmp = OpticalAmplifier('ConstantGain', 20, 5, wavelength); % Can't be ConstantOutputPower here
BWopt = sim.fs;
Npol =2;

PrecUnamp = zeros(size(Lkm));
PrecAPDopt = zeros(size(Lkm));
PrecAmp = zeros(size(Lkm));
OSNRdB = zeros(size(Lkm));
for k = 1:length(Lkm)
    fprintf('L = %d\n', Lkm(k));
    
    Fiber.L = 1e3*Lkm(k);
    
    Hch = Hdac.*Hmod.*Fiber.Himdd(f, wavelength, alpha, 'small signal').*Hrx;
    
    [~, eq] = equalize(eq, [], Hch, mpam, sim, false);
    
    Htot = Hrx.*eq.Hrx.*eq.Hff(sim.f/(mpam.Rs)); % antialiasing + matched + linear equalizer filters
    noiseBW = 1/2*trapz(f, abs(Htot).^2); % filter noise BW (includes noise enhancement penalty)
    
    noiseSTDUnamp = PD.stdNoise(Hrx.*eq.Hrx, eq.Hff(sim.f/(mpam.Rs)), N0, RIN, sim); % thermal + shot + RIN
    
    noiseSTDAmp = @(Plevel) sqrt(noiseBW*N0 + PD.varShot(Plevel, noiseBW) + varRIN(Plevel, noiseBW)... % thermal + shot + RIN
                    + PD.R^2*(OptAmp.varNoiseDD(Plevel/(OptAmp.Gain), noiseBW, BWopt, Npol))); % sig-spont + spont-spont

    PrecUnamp(k) = mpam.sensitivity(noiseSTDUnamp, sim.rexdB, sim.BERtarget);
    PrecAmp(k) = mpam.sensitivity(noiseSTDAmp, sim.rexdB, sim.BERtarget);
    OSNRdB(k) = 10*log10(dBm2Watt(PrecAmp(k))/(2*OptAmp.Ssp*OptAmp.BWref));                
                
    %% APD-based system: shot and thermal noise comparable
    % Gain optimization
    GainAPD = 1:30;
    PrecAPDGain = zeros(size(GainAPD));
    for kk = 1:length(GainAPD)
        APD.Gain = GainAPD(kk);
        Hapd = APD.H(f); % photodiode frequency response
        Hch = Hdac.*Hmod.*Fiber.Himdd(f, wavelength, alpha, 'small signal').*Hapd.*Hrx;  % Channel frequency response
        [~, eq] = equalize(eq, [], Hch, mpam, sim, false);
        
        Htot = Hrx.*eq.Hrx.*eq.Hff(sim.f/(mpam.Rs)); % antialiasing + matched + linear equalizer filters
        noiseBW = 1/2*trapz(f, abs(Htot).^2); % filter noise BW (includes noise enhancement penalty)
        
        noiseSTDAPd = APD.stdNoise(Hrx.*eq.Hrx, eq.Hff(sim.f/(mpam.Rs)), N0, RIN, sim); % thermal + shot + RIN         

        PrecAPDGain(kk) = mpam.sensitivity(noiseSTDAPd, sim.rexdB, sim.BERtarget);
        
        PrecAPDGain(kk) = PrecAPDGain(kk) - 10*log10(APD.Geff); % refer to input power
    end
                
    % Calculate performance at optimal gain
    g = linspace(GainAPD(1), GainAPD(end));
    Prec = interp1(GainAPD, PrecAPDGain, g);    
    [PrecAPDopt(k), idx] = min(Prec);
    APD.Gain = g(idx);
    
    figure(100), hold on, box on
    plot(g, Prec, 'LineWidth', 2)
    plot(APD.Gain, PrecAPDopt(k), '*r')
    xlabel('APD gain')
    ylabel('Receiver sensitivity (dBm)')
    drawnow
    
    Prec = dBm2Watt(PrecAPDopt(k));
    varS = APD.varShot(Prec, noiseBW);
    disp('-- APD')
    fprintf('Optimal gain = %.2f\n', APD.Gain)
    fprintf('Thermal/Shot = %f dB\nRIN/Shot = %f dB\n',...
    10*log10(N0*noiseBW/varS), 10*log10(varRIN(APD.Geff*Prec, noiseBW)/varS))
   
end

PrecAmp = PrecAmp - OptAmp.GaindB; % Refer required power to before the optical amplifier
D = 1e6*abs(Fiber.D(wavelength)*Lkm);

figure(1), hold on, box on
hplot(1) = plot(D, PrecUnamp, '-', 'LineWidth', 2, 'DisplayName', '4-PAM unamplified');
hplot(2) = plot(D, PrecAmp, '--', 'Color', get(hplot(1), 'Color'), 'LineWidth', 2, 'DisplayName', '4-PAM amplified');
hplot(3) = plot(D, PrecAPDopt, ':', 'Color', get(hplot(1), 'Color'), 'LineWidth', 2, 'DisplayName', '4-PAM APD');
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