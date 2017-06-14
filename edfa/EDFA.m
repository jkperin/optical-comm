classdef EDFA
    %% Erbium doped fiber amplification section 
    % References: 
    % [1] C.R. Giles, E. Desurvire, Modeling erbium-doped fiber amplifiers, J. Lightw. Technol. 9 (1991) 271–283. doi:10.1109/50.65886.
       
    properties
        Pump % array of structs containing pump parameters: (P, wavelength, direction, alpha, gain_coeff)
        Signal % array of structs containing pump parameters: (P, wavelength, direction, alpha, gain_coeff)
        L % section length in m
        sat_param=1.5e15 % saturation parameter in m^-1
        BWref=125e9 % reference bandwidth to measure ASE in Hz
    end
    
    properties (Constant, Hidden)
        h = 6.62606957e-34; % Planck
        q = 1.60217657e-19; % electron charge
        c = 299792458;      % speed of light
    end
    
    properties (Constant, Hidden)
        maxL=100; % maximum fiber length used in optimization
    end
    
    methods
        
        function obj = EDFA(Pump, Signal, L, xi)
            %% Constructor
            obj.Pump = Pump;
            obj.Signal = Signal;
            obj.L = L;
            
            if exist('xi', 'var')
                obj.xi = xi;
            end
        end
        
        function [Gain, Psignal_out, Ppump_out] = analytical_gain(self)
            %% Calculate gain by solving analytical equations (24) and (25) of [1]
            Qpump = self.Pump.P./(self.h*self.c./self.Pump.wavelength);
            Qsignal = self.Signal.P./(self.h*self.c./self.Signal.wavelength);
            
            Qin_k = [Qpump Qsignal];
            Qin = sum(Qin_k);
            alpha = [self.Pump.alpha self.Signal.alpha];
            g = [self.Pump.gain_coeff self.Signal.gain_coeff];
            
            % Solve implicit equation (25) of [1] for the total output flux Qout
            Qout0 = Qin;
            [Qout, ~, exitflag] = fzero(@(Qout) Qout - sum(Qin_k.*exp((alpha + g)*(Qin - Qout)/self.sat_param - alpha*self.L)), Qout0);
            
            if exitflag ~= 1
                warning('EDFA/analytical_gain: could not solve for Qout. Simulation exited with exitflag = %d\n', exitflag)
            end
            
            % Calculate the output flux for each individual signal & pump
            % using equation (24) of [1]
            Qout_k = Qin_k.*exp((alpha + g)/self.sat_param*(Qin - Qout) - alpha*self.L);
            
            Ppump_out = Qout_k(1:length(self.Pump.wavelength))*(self.h*self.c./self.Pump.wavelength);
            Psignal_out = Qout_k(length(self.Pump.wavelength)+1:end).*(self.h*self.c./self.Signal.wavelength);
            
            Gain = Psignal_out./self.Signal.P;
        end 
        
        function [ASE, sol] = amplified_power(self)
            [~, Psignal_out, Ppump_out] = self.analytical_gain();
            
            dz = 1;
            z = 0:dz:self.L;
            
            BW = [1520 1570];
            dBW = 1;
            
            nu = 1e-9*(BW(1):dBW:BW(2));
            ASE.wavelength = [nu nu];
            ASE.P = zeros(1, 2*length(nu));
            ASE.direction = ones(size(ASE.P));
            ASE.direction(1:length(nu)) = -1;
            ASE.alpha = 10^(2.6/10)*ones(size(ASE.P));
            ASE.gain_coeff = 10^(3.6/10)*ones(size(ASE.P));
            
            solinit = bvpinit(z, [Ppump_out, Psignal_out, ASE.P]);
            sol = bvp4c(@(z, P) odefun(z, P, self, ASE), @(P0, PL) bcfun(P0, PL, self, ASE), solinit);
            
            Nsignals = length(self.Signal.wavelength);
            Npump = length(self.Pump.wavelength);
            ASE.PoutdBm_b = 10*log10(sol.y(Nsignals + Npump + (1:length(nu)), :)/1e-3);
            ASE.PoutdBm_f = 10*log10(sol.y((Nsignals+Npump+length(nu)+1):end, :)/1e-3);
                        
            function dP = odefun(~, P, Amp, ASE)
                %% Build differential equations dP/dz = odefun(z, P)
                m = 2; % supported number of modes 
                u = [Amp.Pump.direction Amp.Signal.direction ASE.direction].';
                alpha = [Amp.Pump.alpha Amp.Signal.alpha ASE.alpha].';
                g = [Amp.Pump.gain_coeff Amp.Signal.gain_coeff ASE.gain_coeff].';
                lamb = [Amp.Pump.wavelength Amp.Signal.wavelength ASE.wavelength].';
                l = zeros(size(alpha)); % excess loss
                
                n2nt = sum(P.*alpha.*lamb)/(Amp.sat_param*Amp.h*Amp.c + sum(P.*(alpha + g).*lamb));
                
                dP = u.*(alpha + g).*n2nt.*P... % medium gain
                    + u.*g.*n2nt.*m*Amp.h*Amp.c./lamb*Amp.BWref... % ASE
                    - u.*(alpha + l).*P;
            end
            
            function res = bcfun(ya, ~, Amp, ASE)
                %% Calculate residual at boundary conditions
                P = [Amp.Pump.P Amp.Signal.P ASE.P].';
                res = ya-P; % Make Ppump(z = 0) = Pump.P and Psignal(z = 0) = Signal.P
            end
        end
        
        function Lopt = optimal_length(self)
            %% Optimize amplifier length for particular pump/signal configuration.
            % This function uses the analytical_gain() to obtain the gain
            % Length is set to maximize the gain at the first signal
            % channel
            Temp = self;
            [Lopt, ~, exitflag] = fminbnd(@(L) objective(Temp, L), 0, self.maxL);
           
            if exitflag ~= 1
                warning('EDFA/optimal_length: could not solve for Qout. Simulation exited with exitflag = %d\n', exitflag)
            end
            
            function y = objective(Temp, L)
                Temp.L = L;
                Gain = Temp.analytical_gain(); 
                y = -10*log10(Gain(1)); 
            end
        end
        
        function [Psat_pump, Psat_signal] = saturation_power(self)
            %% Calculates the saturation power using the definition of the saturation parameter
            Psat_signal = self.c./self.Signal.wavelength*self.h*self.sat_param./(self.Signal.alpha + self.Signal.gain_coeff);
            Psat_pump = self.c./self.Pump.wavelength*self.h*self.sat_param./(self.Pump.alpha + self.Pump.gain_coeff);
        end
        
    end
end