function Eout = mzm(Ein, Vin, Mod)
%% Dual-polarization I/Q Mach-Zender modulator 
% Inputs:
% - Ein : input electric field
% - Vin : driving signal normalized by Vpi 
% - tx : contains modulators parameters
%   - Mod.Hel : electric frequency response
% Output : Eout is the modulated electric field

% Breaks data into in-phase, quadrature in two pols
Vxi = real(Vin(1, :));
Vxq = imag(Vin(1, :));
Vyi = real(Vin(2, :));
Vyq = imag(Vin(2, :));

% Filter drive waveforms for modulators
Hmod = ifftshift(Mod.Hel);
Vxi = real(ifft(fft(Vxi).*Hmod));
Vxq = real(ifft(fft(Vxq).*Hmod)); 
Vyi = real(ifft(fft(Vyi).*Hmod)); 
Vyq = real(ifft(fft(Vyq).*Hmod));

% Modulate signal fields (each has unit average power)
Vout   = [sin(pi*Vxi/2) + 1i*sin(pi*Vxq/2);...        % x polarization
          sin(pi*Vyi/2) + 1i*sin(pi*Vyq/2)];         % y polarization

Eout =  [Ein/sqrt(2).*Vout(1, :); Ein/sqrt(2).*Vout(2, :)];  % polarization multiplexed signal    