classdef fiber
    % Single-mode fiber
    
    properties
        L % fiber length (m)
        att % attenuation function i.e., alpha = att(lambda) (dB/km)
        D % dispersion function i.e., Dispersion = D(lambda)
    end
    
    properties(Constant)
        S0 = 0.092*1e3;     % dispersion slope (in s/m^3)
        lamb0 = 1310e-9;   % zero-dispersion wavelength
    end
    
    properties(Constant, GetAccess=private)
        c = 299792458;  % speed of light
    end
    
    methods
        function obj = fiber(L, att, D)
            %% Class constructor 
            % Inputs: 
            % - L = fiber length (m)
            % - att (optional, default is no attenuation) function handle of 
            % fiber attenuation (dB/km) as a function of the wavelength. 
            % - D (optional, default is dispersion of standard SMF) = 
            % function handle of fiber chromatic dispersion (s/m^2) as a 
            % function of the wavelength.
            if nargin >=1
                obj.L = L;
            else
                obj.L = 0;
            end
            
            if nargin >= 2
                obj.att = att;
            else % assumes constant attenuation of 0 dB/km
                obj.att = @(lamb) 0;
            end
            
            if nargin == 3
                obj.D = D;
            else % assume SMF28 with zero-dispersion wavelength = 1310nm and slope S0 = 0.092         
                obj.D = @(lamb) fiber.S0/4*(lamb - fiber.lamb0^4./(lamb.^3)); % Dispersion curve
            end           
        end
        
        function link_att = link_attenuation(this, lamb)
            %% Calculate link attenuation in linear units
            % Input:
            % - lamb = wavelength (m)
            link_att = 10^(-this.att(lamb)*this.L/1e4);
        end

        function [Eout, Pout] = linear_propagation(this, Ein, f, lambda)
            %% Linear propagation
            % Perform linear propagation including chromatic dispersion,
            % and attenuation
            % Inputs: 
            % - Ein = input electric field
            % - f = frequency vector (Hz)
            % - lambda = wavelength (m)
            % Outputs: 
            % - Eout, Pout = output electric field and optical power
            % respectively.
            
            if this.L == 0
                Eout = Ein;
                Pout = abs(Ein).^2;
                return
            end
            
            Hele = this.Hdisp(f, lambda);

            Eout = ifft(ifftshift(Hele).*fft(Ein));

            % Received power 
            Eout = Eout*sqrt(this.link_attenuation(lambda));
            
            Pout = abs(Eout).^2;
        end
        
        function Hele = Hdisp(this, f, lambda)
            %% Dispersion frequency response Hele = Eout(f)/Ein(f)
            % Inputs: 
            % - f = frequency vector (Hz)
            % - lambda = wavelength (m)
            % Outputs:
            % - Hele = Eout(f)/Ein(f)
            
            Dispersion =  this.D(lambda);
            beta2 = this.D2beta2(Dispersion, lambda);
            w = 2*pi*f;
            Dw = -1j*1/2*beta2*(w.^2);
            Hele = exp(this.L*Dw);
        end
               
        function Hatt = Hfiber(this, f, tx)
            %% Fiber small-signal frequency response assuming transient chirp dominant
            % This transfer function is for optical power not electic field
            % i.e., Hfiber = Pout/Pin. Moreover, it includes attenuation
            % Inputs:
            % - f = frequency vector (Hz)
            % - tx = transmitter struct. Required fields: lamb (wavelenth
            % in m), and alpha (optional, default zero) (chirp paramter). 
            
            
            beta2 = this.D2beta2(this.D(tx.lamb), tx.lamb);

            % CD frequency response
            theta = -1/2*beta2*(2*pi*f).^2*this.L; % theta = -1/2*beta2*w.^2*L

            if isfield(tx, 'alpha') % if chirp parameter is defined
                alpha = tx.alpha;
            else
                alpha = 0;
            end
            
            H = cos(theta) - alpha*sin(theta);  % fiber small-signal frequency response
            
            % Include attenuation
            Hatt = H.*this.link_attenuation(tx.lamb);
        end  
    end

    methods(Access=private)
        function beta2 = D2beta2(this, D, lamb)
            % Calculates beta2 from dispersion D and wavelength lamb
            beta2 = -D*lamb^2/(2*pi*this.c); 
        end
    end
        

    
end
        