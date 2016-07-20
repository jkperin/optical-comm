function Mod = mzm_frequency_response(ratio, L, f, verbose)
%% Mach-Zehnder frequency response limited by velocity mismatch and loss
% K. P. Ho, Phase-Modulated Optical Communication Systems. New York: Springer, 2005.
% Inputs:
% - ratio: velocity mismatch ratio e.g., 0.98
% - L: Iteractive length in m
% - f: frequency vector
% - verbose (optional, default = false): whether to plot figure
% Output:
% - Mod: struct containing MZM parameters

omega = 2*pi*f;
c = 299792458;      % speed of light
Mod.ratio      = ratio;                                                     % velocity mismatch ratio         
n_r        = 2.15;                                                     % refractive index of coplanar-waveguide (CPW) for TM input light
n_m        = n_r*ratio;                                                % 
Mod.d_12   = (n_m -n_r)/c;                                             % Velocity mismatch difference between the optical and electrical waveguides
a      = 0.01*100*sqrt(abs(omega)/2/pi*1e-9);                      % Microwave attenuation coefficient (electrode loss)
Mod.L      = L;                                                     % Interaction length 
Mod.H    = (1-exp(Mod.L*(-a + 1j*Mod.d_12*omega)))./...          % Freq. response of the optical modulator limited by elec. loss and velocity mismatch
                (a-1j*Mod.d_12*omega)/Mod.L;
Mod.H(isnan(Mod.H)) = 1;                                           % Mod.Hel(f=0) is NaN 
Mod.BW = interp1(abs(Mod.H(f > 0)).^2, f(f > 0), 0.5); 

fprintf('MZM BW = %.2f GHz\n', Mod.BW/1e9);

Mod.H = Hgrpdelay(Mod.H, f); % remove group delay

if exist('verbose', 'var') && verbose
    figure(111), hold on
    plot(f/1e9, abs(Mod.H).^2)
    plot(f([1 end])/1e9, [0.5 0.5], 'k')
    plot(Mod.BW/1e9, 0.5, 'xk', 'MarkerSize', 6)
    xlabel('Frequency (GHz)')
    ylabel('|H(f)|^2')
    legend('MZM frequency response', sprintf('BW = %.2f GHz', Mod.BW/1e9));
    title('MZM frequency response')
end
    
    