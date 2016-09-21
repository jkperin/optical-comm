clear, clc, close all

c = 299792458;  % speed of light
Fspacing = 100e9; % Hz
L = 40; % km
receiverSensitivity = -34; % dBm
maxD = 9e3*6e-6; % 9 km at 1250e-6 leads to 5-dB penalty
att = 0.35; % dB/km
loss = 7;
Pmax = 9.8; % dBm

lamb0 = 1310e-9;
f0 = c/lamb0;
Dlamb = c/f0 - c/(f0 + Fspacing)

%% 100 GHz spacing
disp('> 100-GHz spacing')
Nch = 10;
Disp = 0;
while all(abs(Disp)*L*1e3 <= maxD) 
    [lambchs, Disp] = optimal_wdm_channels(Dlamb, Nch);
    Nch = Nch + 1;
end
Nch  = Nch - 1

Plaunch = Pmax - 10*log10(Nch)
Margin = Plaunch - loss - att*L - receiverSensitivity

disp('> 200-GHz spacing')
Nch = 10;
Disp = 0;
while all(abs(Disp)*L*1e3 <= maxD) 
    [lambchs, Disp] = optimal_wdm_channels(2*Dlamb, Nch);
    Nch = Nch + 1;
end
Nch  = Nch - 1

Plaunch = Pmax - 10*log10(Nch)
Margin = Plaunch - loss - att*L - receiverSensitivity





