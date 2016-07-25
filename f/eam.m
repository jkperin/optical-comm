function [Et, Pt] = eam(Ein, Vin, Mod, f)
%% Electro absorption modulator including the following effects
% 1. frequency response
% 2. Modulator nonlinearity (not implemented)
% 3. Adds chirp if Mod.alpha is defined

% Apply frequency response of the modulator
Pt = abs(Ein).^2./mean(abs(Ein).^2).*real(ifft(fft(Vin).*ifftshift(Mod.H(f))));

% Modulator Nonlinearity
%## To be implemented
    
% Enforce modulator extinction ratio
% if isfield(Mod, 'rexdB')
%     rex = 10^(-abs(Mod.rexdB)/10);
%     Padd = (rex*max(Pt) - min(Pt))/(1 - rex);
%     Pt = Pt + Padd;
% end

Pt(Pt <= 0) = 0; 

assert(all(Pt >= 0), 'eam: Negative power')

Et = sqrt(Pt);

if isfield(Mod, 'alpha') && Mod.alpha ~= 0
    %% Add transient chirp
    % Calculate electric field including chirp
    dphi = Mod.alpha/2*log(Pt + realmin);         % Phase variation due to chirp (only transient chirp is considered)
    % realmin = smallest positive normalized floating-point number.
    % Just to avoid log(0)

    Et = Et.*exp(1j*dphi);
end

