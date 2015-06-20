classdef fiber
    properties
        L % fiber length (m)
        att % attenuation function i.e., alpha = att(lambda) (dB/km)
        D % dispersion function i.e., Dispersion = D(lambda)
    end
    
    properties(Constant, GetAccess=private)
        c = 299792458;  % speed of light
    end
    
    methods
        %% Constructor
        function obj = fiber(L, att, D)
            obj.L = L;
            
            if nargin >= 2
                obj.att = att;
            else % assumes constant attenuation of 0.35 dB/km
                obj.att = @(lamb) 0.35;
            end
            
            if nargin == 3
                obj.D = D;
            else % assume SMF28 with zero-dispersion wavelength = 1310nm and slope S0 = 0.092
                S0 = 0.092*1e3;     % dispersion slope (in s/m^3)
                lamb0 = 1310e-9;   % zero-dispersion wavelength (also used to convert D to beta2)

                obj.D = @(lamb) S0/4*(lamb - lamb0^4./(lamb.^3)); % Dispersion curve
            end           
        end
        
        %% Linear propagation
        function Eout = linear_propagation(this, Ein, f, lambda)
            Dispersion =  fiber.D(lambda);
            beta2 = this.D2beta2(Dispersion, lambda);
            w = 2*pi*ifftshift(f);
            
            % Electric field after fiber propagation
            Dw = -1j*1/2*beta2*(w.^2);
            Eout = ifft(exp(fiber.L*Dw).*fft(Ein));

            % Received power 
            Eout = Eout*sqrt(10^(-fiber.att(lambda)*fiber.L/1e4));
        end       
    end
    
    methods(Access=private)
        function beta2 = D2beta2(this, D, lamb)
            beta2 = -D*lamb^2/(2*pi*this.c); 
        end
    end
        

    
end
        