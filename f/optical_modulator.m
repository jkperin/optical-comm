% 1. frequency response
% 2. Modulator nonlinearity (not implemented)
% 3. Adds intensity noise (if sim.RIN = true)
% 4. Adds chirp if tx.alpha is defined

function [Et, Pt] = optical_modulator(xt, tx, sim)
    %% Apply frequency response of the modulator
    if isfield(tx, 'modulator') % if frequency response is defined
        Hmod = tx.modulator.H(sim.f).*exp(1j*2*pi*sim.f*tx.modulator.grpdelay);
        Pt = tx.kappa*real(ifft(ifftshift(Hmod).*fft(xt)));  % Note: real() is used to remove residual imag
                                                    % part that appears due to numerical error   
    else
        Pt = tx.kappa*xt;
    end

    %% Modulator Nonlinearity
    %## To be implemented
 
    %% Add intensity noise
    if isfield(sim, 'RIN') && sim.RIN
        % Calculate RIN PSD, which depends on the instantaneous power
        Srin = 10^(tx.RIN/10)*Pt.^2;

        % noise. Srin is two-sided psd
        wrin = sqrt(Srin.*sim.fs).*randn(size(Pt));
        
        Pt = Pt + wrin;
    end 
    
    % Ensures signal is non-negative
    Pt(Pt < 0) = 0;
    
    if isfield(tx, 'alpha')
        %% Add chirp
        % Calculate electric field including chirp
        dphi = tx.alpha/2*log(Pt + realmin);         % Phase variation due to chirp (only transient chirp is considered)
        % realmin = smallest positive normalized floating-point number.
        % Just to avoid log(0)
        
        Et = sqrt(Pt).*exp(1j*dphi);
    else
        Et = sqrt(Pt);
    end       
   
end