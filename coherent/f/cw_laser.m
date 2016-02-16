function Eout = cw_laser(laser, sim)
%% Generates CW light source in 1 pol. Includes intensity noise, phase noise, and frequency shift from reference frequency
% Pout = output power in W
% laser = contains laser properties
% - PdBm : output power in dBm
% - RIN in dB/Hz
% - linewidth in Hz
% - freqshift (optional) in Hz : frequency shift from reference frequency
% sim = simulation parameters

Pout = dBm2Watt(laser.PdBm)*ones(1, sim.N);

%% ============ Intensity Noise ==========
if isfield(sim, 'RIN') && sim.RIN
    if isfield(laser, 'RIN_shapedB') % RIN noise modeled by PSD   
        if length(Pout) == length(laser.RIN_shapedB)
            RIN_shape = 10.^(laser.RIN_shapedB/10);
        else
            RIN_shape = 0;
        end

        wrin = sqrt(Pout.^2*sim.fs).*randn(size(sim.f));

        wrin = real(ifft(fft(wrin).*ifftshift(sqrt(RIN_shape)))); 
    else % white noise approximation
        % Calculate RIN PSD, which depends on the instantaneous power
        Srin = 10^(laser.RIN/10)*Pout.^2; 

        % Intensity noise. Srin is one-sided psd
        wrin = sqrt(Srin.*sim.fs/2).*randn(size(Pout));
    end
    
    Eout = sqrt(Pout + wrin); % adds intensity noise;
else
    Eout = sqrt(Pout);
end

%% ============ Phase Noise ==========
if isfield(sim, 'phase_noise') && sim.phase_noise                                                                        
    pn_std = sqrt(2*pi*laser.linewidth/sim.fs);                               % sqrt of variance of phase difference over one symbol
    % Generate random received phase following a Wiener process ---------------
    initial_phase    = pi*(2*rand(1)-1); % [-pi, pi]
    dtheta      = [0 pn_std*randn(1, sim.N-1)];                                  % i.i.d.Gaussian random variables with zero mean and variance sigma_p^2
    phase_noise = initial_phase + cumsum(dtheta, 2);
    
    Eout = Eout.*exp(1j*phase_noise); % adds phase noise
end

% If freqshift is defined adds frequency shift
if isfield(laser, 'freqshift')
    Eout = freqshift(Eout, sim.t, laser.freqshift);
end
