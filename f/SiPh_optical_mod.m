function Eout = SiPh_optical_mod(Ein, Vin, Mod)
%% Silicon-Photonics dual polarization dual quadrature optical modulator 
% Inputs:
% - Ein : input electric field
% - Vin : driving signal
% - Mod : contains modulators parameters
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
Vout   = 2^(-1/4)*[Vxi + 1i*Vxq;...        % x polarization
                    Vyi + 1i*Vyq];         % y polarization

Eout =  1/sqrt(2)*[Ein.*Vout(1, :); Ein.*Vout(2, :)];  % polarization multiplexed signal    