classdef fiber < handle
    %% Single-mode fiber
    properties
        L % fiber length (m)
        att % attenuation function i.e., alpha = att(lambda) (dB/km)
        D % dispersion function i.e., Dispersion = D(lambda)
        PMD = false; % whether PMD is included in simulations
        gamma = 1.4e-3 % nonlinear coefficient in 1/(W*m)
        meanDGDps = 0.1; % mean DGD  (ps/sqrt(km))
        PMD_section_length = 1e3  % Section length for simulating PMD (m)
        PDL = 0 % polarization dependent loss (dB). Here, it indicates how much the y pol will be attenuated with respect to x pol
    end
    
    properties(Dependent)
        tauDGD % total differential group delay (s)
    end
        
    properties(GetAccess=protected)
        JonesMatrix % Jones Matrix 
    end
    
    properties(Constant) % These parameters are only used for the default value of D(lamb)
        S0 = 0.092*1e3;    % dispersion slope (in s/m^3)
        lamb0 = 1310e-9;   % zero-dispersion wavelength
    end
    
    properties(Constant, Hidden)
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
        
        function Fibertable = summary(self, lamb)
            %% Generate table summarizing class values
            disp('-- Fiber class parameters summary:')
            lambnm = lamb*1e9;
            rows = {'Length'; sprintf('Attenuation at %.2f nm', lambnm);...
                sprintf('Total dispersion at %.2f nm', lambnm); 'PMD included?';...
                'Total DGD'; 'Polarization dependent loss'};
            Variables = {'L'; 'att'; 'DL'; 'PMD'; 'tauDGD'; 'PDL'};
            Values = [self.L/1e3; self.att(lamb); self.D(lamb)*self.L*1e3;...
                self.PMD; self.tauDGD*1e12; self.PDL];
            Units = {'km'; 'dB/km'; 'ps/nm'; ''; 'ps'; 'dB'};

            Fibertable = table(Variables, Values, Units, 'RowNames', rows);
        end
             
        %% Get methods
        function tauDGD = get.tauDGD(self)
            tauDGD = double(self.PMD)*self.meanDGDps*1e-12*sqrt(self.L/1e3); % corresponds to total 
        end
        
        %% Main Methods
        function [Ncd, Npmd] = Ntaps(self, Rs, ros, lambda)
            %% Estimated number of taps in DSP required to compensate for CD and PMD in coherent detection link
            % Based on Ip, E., & Kahn, J. M. (2007). Digital equalization of chromatic dispersion and polarization mode dispersion. 
            % Journal of Lightwave Technology, 25(8), 2033–2043.
            Ncd = 2*pi*abs(self.beta2(lambda))*self.L*Rs^2*ros; % eq (35)
            Npmd = self.tauDGD*ros*Rs; % eq (36)
        end
        
        function b2 = beta2(this, lamb)
            %% Calculates beta2 at wavelength lamb
            b2 = -this.D(lamb).*(lamb.^2)/(2*pi*this.c); 
        end    
        
        function Leff = effective_length(self, lamb)
            alpha = log(10)*self.att(lamb)/1e4; % 1/m
            Leff = (1 - exp(-alpha*self.L))/alpha;
        end
        
        function [link_att, link_attdB] = link_attenuation(this, lamb)
            %% Calculate link attenuation in linear units
            % Input:
            % - lamb = wavelength (m)
            link_att = 10.^(-this.att(lamb)*this.L/1e4); % exp(-alpha*L)
            link_attdB = this.L/1e3*this.att(lamb);
        end
        
        function Eout = linear_propagation(this, Ein, f, lambda)
            %% Linear propagation including chromatic dispersion and first-order PMD, if PMD flag is set
            % Perform linear propagation including chromatic dispersion,
            % and attenuation
            % Inputs: 
            % - Ein = input electric field
            % - f = frequency vector (Hz)
            % - lambda = wavelength (m)
            % Outputs: 
            % - Eout, Pout = output electric field and optical power
            % respectively.
            
            % back-to-back case
            if this.L*this.D(lambda) == 0
                Eout = Ein;
                return
            end
            
            % Make sure that vectors are in the form [1 x N] for 1 pol or [2 x N] for two pols
            Ein = this.enforceDimConvention(Ein);
            f = this.enforceDimConvention(f);
         
            two_pols = (size(Ein, 1) == 2);

            % Frequency domain
            Einf = fftshift(fft(Ein, [], 2), 2);
            
            % PMD
            if this.PMD
                if isempty(this.JonesMatrix) % only calculates Jones Matrix if it doesn't already exist
                    this.generateJonesMatrix(2*pi*f); % 2 x 2 x length(f)
                end
                
                if not(two_pols)
                    Einf = [Einf; zeros(size(Einf))];
                    two_pols = true;
                end

                for k = 1:length(f)
                    Einf(:, k) = this.JonesMatrix(:,:,k)*Einf(:, k);
                end
            end

            % Chromatic dispersion
            Hele = this.Hdisp(f, lambda);
            if two_pols
                Eout = Einf;
                Eout(1, :) = ifft(ifftshift(Hele.*Einf(1, :)));
                Eout(2, :) = ifft(ifftshift(Hele.*Einf(2, :)));
            else
                Eout = ifft(ifftshift(Hele.*Einf));
            end

            % PDL
            if two_pols && this.PDL ~= 0
                a = 10^(-this.PDL/10);
                Eout(2, :) = a*Eout(2, :);
            end
            
            % Attenuation
            Eout = Eout*sqrt(this.link_attenuation(lambda));
        end
        
        function Hele = Hdisp(this, f, lambda)
            %% Dispersion frequency response Hele(f) = Eout(f)/Ein(f)
            % Inputs: 
            % - f = frequency vector (Hz)
            % - lambda = wavelength (m)
            % Outputs:
            % - Hele = Eout(f)/Ein(f)
            
            beta2 = this.beta2(lambda);
            w = 2*pi*f;
            Dw = -1j*1/2*beta2*(w.^2);
            Hele = exp(this.L*Dw);
        end
        
        function h = hdisp(self, t, lambda)
            %% Dispersion impulse response inverseFourier(Hdisp(f))
            % Inputs: 
            % - f = frequency vector (Hz)
            % - lambda = wavelength (m)

            b = 1/2*self.beta2(lambda)*self.L;
            h = sqrt(-pi*1j/b)/(2*pi)*exp(1j*t.^2/(4*b));
            % h = sqrt(1/(4*pi*b))*(cos(t.^2/(4*b) - pi/4) + 1j*sin(t.^2/(4*b) - pi/4));
        end
               
        function Hf = Himdd(self, f, wavelength, alpha, type)
            %% Fiber frequency response for an IM-DD system with transient chirp dominant
            % This transfer function is for optical power not electic field
            % i.e., Hfiber(f) = Pout(f)/Pin(f).
            % Inputs:
            % - f: frequency vector (Hz)
            % - wavelength: wavelength (m)
            % - alpha (optional, default alpha = 0): chirp parameter with
            % sign convention such that for DML alpha > 0
            % - type (optional, default type = 'small signal')
            % in m), and alpha (optional, default zero) (chirp paramter). 
 
            if not(exist('alpha', 'var')) % if chirp parameter is not defined
                alpha = 0;
            end
            
            beta2 = self.beta2(wavelength);
            theta = -1/2*beta2*(2*pi*f).^2*self.L; % theta = -1/2*beta2*w.^2*L
            
            if exist('type', 'var') && strcmpi(type, 'large signal')
                %% Large signal
                mIM = 0.7; % modulation index is set to 70%
                Dphi = pi/2; % i.e., transient chirp dominant
                mFM = alpha/2*mIM; 
                u = 2*mFM*sin(theta);
                Hf = cos(theta).*(besselj(0, u) - besselj(2, u)*exp(1j*Dphi)) - 2*exp(1j*Dphi)/(1j*mIM)*besselj(1, u);  % fiber large-signal frequency response                 
            elseif not(exist('type', 'var')) || strcmpi(type, 'small signal')
                %% Small signal
                Hf = cos(theta) - alpha*sin(theta);  % fiber small-signal frequency response
            else
                error('fiber/Hf: undefined type of fiber frequency response')
            end
        end
        
        function tau = calcDGD(self, omega)
            %% Calculate differential group delay from Jones Matrix
            if ~self.PMD 
                tau = zeros(size(omega));
                warning('fiber/calcDGD: PMD is disable')
                return
            end
                
            if isempty(self.JonesMatrix)
                tau = [];
                warning('fiber/calcDGD: Jones Matrix was not calculated yet')
                return
            end
            
            tau = zeros(1,length(omega));
            dw = abs(omega(1)-omega(2));
            for m = 1:length(omega)-1;
                tau(m) = 2/dw*sqrt(det(self.JonesMatrix(:,:,m+1)-self.JonesMatrix(:,:,m)));
            end
            tau(end) = tau(end-1);
         end
        
    end  

    methods(Access=private)
        function V = enforceDimConvention(~, V)
            %% Ensures that V is in the form [1 x N] for 1 pol or [2 x N] for 2 pols
            if size(V, 1) > size(V, 2)
                V = V.';
            end
        end
        
        function M = generateJonesMatrix(self, omega)
            %% Function to generate Jones Matrix, modified from Milad Sharif's code
            Nsect = ceil(self.L/self.PMD_section_length);
                       
            dtau = self.tauDGD/sqrt(Nsect);

            M = zeros(2,2,length(omega));
            M(:, :, 1) = randomRotationMatrix(); % rotation to an arbritary polarizataion state
            
            for k = 1:Nsect
                U = randomRotationMatrix();

                for m = 2:length(omega)
                    Dw = [exp(1j*dtau*omega(m)/2), 0; 0, exp(-1j*dtau*omega(m)/2)]; % Birefringence matrix
                    M(:,:,m) = U*Dw*U';
                end
            end
            
            self.JonesMatrix = M;
            
            function U = randomRotationMatrix()
                phi = rand(1, 3)*2*pi;
                U1 = [exp(-1j*phi(1)/2), 0; 0 exp(1j*phi(1)/2)];
                U2 = [cos(phi(2)/2) -1j*sin(phi(2)/2); -1j*sin(phi(2)/2) cos(phi(2)/2)];
                U3 = [cos(phi(3)/2) -sin(phi(3)/2); sin(phi(3)/2) cos(phi(3)/2)];

                U = U1*U2*U3;
%                U = [cos(phi), -sin(phi); sin(phi), cos(phi)];
            end            
        end
    end    
end
