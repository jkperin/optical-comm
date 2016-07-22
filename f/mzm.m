function Eout = mzm(Ein, Vin, Mod)
%% Mach-Zehnder modulator
% - Sual-pol I/Q modulator if Vin is 2 x N and complex
% - Single-pol I/Q modulator if Vin is 1 x N and complex
% - Intensity modulator if Vin is 1 x N and real
% Inputs:
% - Ein : input electric field
% - Vin : driving signal normalized by Vpi/2. Input signal must already
% have proper biasing and voltage swing
% - tx : contains modulators parameters
%   - Mod.Hel : electric frequency response
% Output : Eout is the modulated electric field
Hmod = ifftshift(Mod.H);    

if size(Vin, 1) == 2 && not(isreal(Vin)) 
%% Dual-pol I/Q
    % Breaks data into in-phase, quadrature in two pols
    Vxi = real(Vin(1, :));
    Vxq = imag(Vin(1, :));
    Vyi = real(Vin(2, :));
    Vyq = imag(Vin(2, :));

    % Filter driving signals
    Vxi = real(ifft(fft(Vxi).*Hmod));
    Vxq = real(ifft(fft(Vxq).*Hmod)); 
    Vyi = real(ifft(fft(Vyi).*Hmod)); 
    Vyq = real(ifft(fft(Vyq).*Hmod));

    % Modulate signal fields (each has unit average power)
    Enorm = 1/sqrt(sqrt(2)); % normalize so that E(|Vout|^2) = 1
    Vout   = Enorm*[sin(pi*Vxi/2) + 1i*sin(pi*Vxq/2);...        % x polarization
                    sin(pi*Vyi/2) + 1i*sin(pi*Vyq/2)];         % y polarization

    Eout =  [Ein/sqrt(2).*Vout(1, :); Ein/sqrt(2).*Vout(2, :)];  % polarization multiplexed signal    
elseif size(Vin, 1) == 1 && not(isreal(Vin)) 
%% Single-pol I/Q
    % Breaks data into in-phase, quadrature
    Vxi = real(Vin(1, :));
    Vxq = imag(Vin(1, :));

    % Filter driving signals
    Vxi = real(ifft(fft(Vxi).*Hmod));
    Vxq = real(ifft(fft(Vxq).*Hmod)); 

    % Modulate signal fields (each has unit average power)
    Enorm = 1/sqrt(sqrt(2)); % normalize so that E(|Vout|^2) = 1
    Vout   = Enorm*(sin(pi*Vxi/2) + 1i*sin(pi*Vxq/2));

    Eout =  Ein.*Vout(1, :);  % polarization multiplexed signal
elseif size(Vin, 1) == 1 && isreal(Vin) 
%% Single-pol intensity modulator
    % Breaks data into in-phase, quadrature in two pols
    Vx = Vin(1, :);

    % Filter driving signal
    Vx = real(ifft(fft(Vx).*Hmod));

    % Modulate signal fields (each has unit average power)
    Enorm = 1; % normalize so that E(|Vout|^2) = 1, if Vx in [-1, 1]
    Vout   = Enorm*sin(pi*Vx/2);

    Eout =  Ein.*Vout(1, :);  % polarization multiplexed signal
else
    error('mzm: Invalid driving signal. Vin must be either 1 x N (real or complex) or 2 x N (complex)')
end