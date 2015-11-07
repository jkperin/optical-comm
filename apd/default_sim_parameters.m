%% SOA and APD default simulation parameters

%% Shot and RIN included
sim.shot = true; % include shot noise
sim.RIN = true; % include RIN noise
sim.BERtarget = 1.8e-4; % target BER

%% Transmitter
tx.lamb0 = 1310e-9; % zero-dispersion wavelength
tx.lamb = 1270e-9; % worst-case wavelength (this may have to change)
tx.alpha = 2; % chirp parameter
tx.RIN = -150;  % dB/Hz
tx.rexdB = -10;  % extinction ratio in dB. Defined as Pmin/Pmax

% Modulator frequency response
tx.modulator.fc = 30e9; % modulator cut off frequency

%% Receiver
rx.N0 = (30e-12).^2; % thermal noise psd
rx.Id = 10e-9; % dark current 
rx.R = 1; % responsivity

%% Fiber
fiber.att = @(lamb) 0.35; % dB/km attenuation constant with wavelength 
fiber.S0 = 0.092*1e3;     % dispersion slope (in s/m^3) (SSMF)
fiber.D = @(lamb) fiber.S0/4*(lamb - fiber.lamb0^4./(lamb.^3)); % Dispersion curve
% assume SMF28 with zero-dispersion wavelength = 1310nm and slope S0 = 0.092   

%% Equalization
% Receiver filter = matched filter
% Equalizer = Fixed time-domain linear equalizer

%% Other effects
% No quantization
% No modulator nonlinearities
