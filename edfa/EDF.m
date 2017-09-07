classdef EDF
    %% Single-mode Erbium-doped fiber
    % EDF parameters are loaded from function edf_selection.m
      
    %% Numerical modeling of EDF physics assumes the standard confined doping (SCD) model [1], [2, Chap. 1], [3, Chap. 1]
    % This model makes the following assumptions
    % 1. Pump, signal, and ASE propagate in the fiber fundamental mode
    % 2. The gain medium is homogeneously broadened
    % 3. ASE is generated in both polarizations
    % 4. The Er-doping is confined to the fiber core, which permits analytical 
    % solution of the transversal intergrals. This last assumption is what
    % differentiates the SCD model from the more general model presented in [2, Chap. 1]
    
    %% References: 
    % [1] C. Giles, E. Desurvire, "Modeling erbium-doped fiber amplifiers,"
    % J. Lightw. Technol. 9 (1991) 271–283.
    % [2] E. Desurvire. "Erbium-Doped Fiber Amplifiers: Principles and Applications," 1994.
    % [3] E. Desurvire, D. Bayart, B. Desthieux, and S. Bigo. "Erbium-Doped 
    % Fiber Amplifiers: Device and System Developments". 2002.
    % [4] P.C. Becker, N.A. Olsson, and J.R. Simpson. "Erbium-doped fiber
    % amplifiers: fundamentals and technology," 1999
    properties
        type='giles_ge:silicate' % which fiber to use. 
        % Fibers available:  
        % - {'giles_ge:silicate' (default), 'giles_al:ge:silicate'}. Extracted from [1, Fig. 2].
        % - {'principles_type1', 'principles_type2', 'principles_type3'}. 3 types of Er3+ in alumino-germanosilicate glass. Figs. 4.20, 4.21, and 4.22 of [2]
        L % fiber length (m)
        excess_loss = 0 % (dB/m) excess loss due to splices for instance. 
        core_radius = 1.4e-6 % Fiber core radius. e.g, 1.2 um in [1, Table 1], 1.4 um in [4, pg. 156]
        doping_radius = 1.2e-6 % Er3+ core radius. e.g., 1.2 um in [1, Table 1], 1.05um in [4, pg 156]
        rho0 = 0.7e19; % Er3+ concentraction (cm^3) [4, pg 156]
        NA = 0.28 % numerical aperture, value taken from [4, pg. 156]
        tau = 10e-3; % metastable lifetime in s        
        Nmode = 2 % number of modes (default = 2 for two polarizations)
    end
    
    properties (Dependent)
        param % physical parameters from fiber (loaded from edf_selection.m)
    end
    
    properties (Constant, Hidden)
        abs_cs980nm = 2.7e-25 % absorption cross section near 980 nm (m^2). 
        % Value obtainend from [4, pg 154]. Cross section curves in edf_select.m are only valid near 1550 nm
        % Emission cross-section is assumed 0, so that two-level system can
        % be used
        maxL = 30; % maximum fiber length. Used to limit simulation
        h = 6.62606957e-34; % Planck
        q = 1.60217657e-19; % electron charge
        c = 299792458;      % speed of light
    end
        
    methods
        function obj = EDF(L, type)
            %% Constructor
            obj.L = L;
            if exist('type', 'var')
                obj.type = type;
            end
        end
        
        %% Get methods
        function param = get.param(self)
            %% Load EDF parameters
            param = edf_selection(self.type);
        end
                       
        %% Semi-analytical models
        function [GaindB, Psignal_out, Ppump_out] = semi_analytical_gain(self, Pump, Signal)
            %% Calculate gain by solving implicit equation (25) of [1]
            % This model assumes
            % 1. Gamma (ovelap integral) is constant
            % 2. amplifier is not saturated by ASE
            % 3. Excess loss is negligible
            % This model does not depend on the direction of propagation
            % between pump and signal
            % Inputs:
            % - Pump: instance of class Channels representing the pump
            % - Signal: instance of class Channels representing the signals           
            Qpump = Pump.P./self.Ephoton(Pump.wavelength);
            Qsignal = Signal.P./self.Ephoton(Signal.wavelength);
            
            Qin_k = [Qpump Qsignal]; % concat pump and signal
            Qin = sum(Qin_k);        % total input power
            
            % Absorption and gain coefficients
            lamb = [Pump.wavelength Signal.wavelength];
            alpha = self.absorption_coeff(lamb); % (1/m)
            g = self.gain_coeff(lamb); % (1/m)
            xi = self.sat_param(lamb);
                       
            % Solve implicit equation (25) of [1] for the total output flux Qout
            Qout0 = Qin; % starting point
            [Qout, ~, exitflag] = fzero(@(Qout) Qout - sum(Qin_k.*exp((alpha + g)*(Qin - Qout)./xi - alpha*self.L)), Qout0);
            
            if exitflag ~= 1
                warning('EDF/analytical_gain: could not solve for Qout. Simulation exited with exitflag = %d\n', exitflag)
            end
            
            % Calculate the output flux for each individual signal & pump
            % using equation (24) of [1]
            Qout_k = Qin_k.*exp((alpha + g).*(Qin - Qout)./xi - alpha*self.L);
            
            Ppump_out = Qout_k(1:Pump.N).*self.Ephoton(Pump.wavelength);
            Psignal_out = Qout_k(Pump.N+1:end).*self.Ephoton(Signal.wavelength);
            
            GaindB = 10*log10(Psignal_out./Signal.P);
        end 
        
        % Pump dominant (incomplete!)
        function [GaindB, Psignal_out, Ppump_out] = gain_approx(self, Pump, Signal)        
            Qin_p = Pump.P./self.Ephoton(Pump.wavelength);
            Qin_k = Signal.P./self.Ephoton(Signal.wavelength);
            
            Qin = Qin_p;        % pump is much larger than signals
            
            % Absorption and gain coefficients
            alphap = self.absorption_coeff(Pump.wavelength);
            gp = self.gain_coeff(Pump.wavelength);
            xip = self.sat_param(Pump.wavelength);
            alpha = self.absorption_coeff(Signal.wavelength); % (1/m)
            g = self.gain_coeff(Signal.wavelength); % (1/m)
            xi = self.sat_param(Signal.wavelength);
                       
            % Solve implicit equation (25) of [1] for the total output flux Qout
            Qout0 = Qin; % starting point
            [Qout, ~, exitflag] = fzero(@(Qout) Qout - Qin_p.*exp((alphap + gp)*(Qin - Qout)./xip - alphap*self.L), Qout0);
            
            if exitflag ~= 1
                warning('EDF/analytical_gain: could not solve for Qout. Simulation exited with exitflag = %d\n', exitflag)
            end
            
            % Calculate the output flux for each individual signal & pump
            % using equation (24) of [1]
            Qout_k = Qin_k.*exp((alpha + g).*(Qin - Qout)./xi - alpha*self.L);
            
            Ppump_out = Qout.*self.Ephoton(Pump.wavelength);
            Psignal_out = Qout_k.*self.Ephoton(Signal.wavelength);
            
            GaindB = 10*log10(Psignal_out./Signal.P);
        end        
        
        
        %% Analytical model
        function GaindB = analytical_gain(self, Pump, Signal)
            %% Calculate gain assuming that fiber is uniformly inverted. Does not depend on pump power
            alphap = self.absorption_coeff(Pump.wavelength);
            gp = self.gain_coeff(Pump.wavelength);
            
            alphas = self.absorption_coeff(Signal.wavelength);
            gs = self.gain_coeff(Signal.wavelength);
            
            % ln(Gain)
            g = self.L*((alphas+gs)*alphap/(alphap+gp) - alphas);
            
            GaindB = 10*log10(exp(g));
        end
           
        function nsp = analytical_excess_noise(self, Pump, Signal)
            %% Analytical expression for excess noise [1, eq. (31)]
            [~, gp] = self.gain_coeff(Pump.wavelength);
            [~, alphap] = self.absorption_coeff(Pump.wavelength);
            
            [~, gs] = self.gain_coeff(Signal.wavelength);
            [~, alphas] = self.absorption_coeff(Signal.wavelength);
            
            nsp = 1./(1 - gp*alphas./(alphap*gs));
        end
        
        function Pase = analytical_ASE_PSD(self, Pump, Signal)
            %% Analytical expression for ASE PSD [1, eq. (32)]
            nsp = analytical_excess_noise(self, Pump, Signal);
            GdB = semi_analytical_gain(self, Pump, Signal);
            Pase = 2*nsp.*(10.^(GdB/10)-1).*self.Ephoton(Signal.wavelength);
        end
        
        function Lopt = optimal_length(self, Pump, Signal)
            %% Optimize EDF length to maximize mean gain for signals with non-zero power
            % This function uses the semi_analytical_gain() to obtain the gain
            Temp = self;
            [Lopt, ~, exitflag] = fminbnd(@(L) objective(Temp, L), 0, self.maxL);
           
            if exitflag ~= 1
                warning('EDFA/optimal_length: could not solve for Qout. Simulation exited with exitflag = %d\n', exitflag)
            end
            
            function y = objective(Temp, L)
                Temp.L = L;
                GaindB = Temp.semi_analytical_gain(Pump, Signal); 
                y = -mean(10.^(GaindB(Signal.P ~= 0)/10)); % objective is to maximize mean gain over channels with non-zero power 
            end
        end
        
        %% Numerical models 
        function [GaindB, Ppump_out, Psignal_out, Pase, sol]...
                = two_level_system(self, Pump, Signal, ASEf, ASEb, BWref, Nsteps)
            %% Calculate Gain and noise PSD by solving coupled nonlinear first-order differential equations that follow from the SCD model
            % This assumes that amplifier is pumped as a two-level system.
            % This is true for pump near wavelength 1480 nm and an
            % approximation for pump near 980 nm. 
                                              
            lamb = [Pump.wavelength, Signal.wavelength, ASEf.wavelength, ASEb.wavelength].'; % wavelength
            u = [Pump.u, Signal.u, ASEf.u, ASEb.u].'; % propation direction
            g = self.gain_coeff(lamb); % Gain coefficient
            alpha = self.absorption_coeff(lamb); % Absorption coefficient
            
            % ASE term is only included in the forward and backward ASE channels
            ASEselect = ones(size(lamb));
            ASEselect(1:(Pump.N+Signal.N)) = 0;
            
            % Photon energy times saturation parameter
            h_nu_xi = self.Ephoton(lamb).*self.sat_param(lamb); % h*nu*xi
            
            % Solver
            options = bvpset('Vectorized', 'on');
            z = linspace(0, self.L, Nsteps).';
            solinit = bvpinit(z, [Pump.P Signal.P ASEf.P ASEb.P].'); % initial guess
            sol = bvp4c(@(z, P) odefun(z, P, self, lamb, h_nu_xi, u, g, alpha, BWref),... % differential equation
                @(P0, PL) bcfun(P0, PL, Pump, Signal, ASEf, ASEb),... % boundary conditions
                solinit, options); % initial guess
            
            if any(sol.y < 0)
                warning('EDF/two_level_system: solution contains negative power')
            end
 
            if strcmpi(Pump.direction, 'forward')
                Ppump_out = sol.y(1:Pump.N, end).';
            else
                Ppump_out = sol.y(1:Pump.N, 1).';
            end
                
            Psignal_out = sol.y(Pump.N + (1:Signal.N), end).';
            Pase = sol.y(Pump.N + Signal.N + (1:Signal.N), end).'; % only forward ASE
            GaindB = 10*log10(Psignal_out./Signal.P);            
            
            function dP = odefun(~, P, edf, lamb, h_nu_xi, u, g, alpha, BWref)
                %% Build differential equations dP/dz = odefun(z, P)
                % Inputs:
                % - z (not used): distance (m)
                % - P: power, including pump, signal, forward ASE, and backward ASE (W)
                % - lamb: wavelengths, including pump, signal, forward ASE, and backward ASE (W)
                % - h_nu_xi: photon energy times saturation parameter
                % - u: propgation direction (+1 if forward, -1 if backward)
                % - g: gain coefficient at wavelengths given in lamb
                % - alpha: absorption coefficient at wavelengths given in lamb
                % - BWref: bandwidth over which to measure ASE
                
                P(P < 0) = 0; % non-negative constraint
                
                n2 = sum(P.*alpha./h_nu_xi)./(1 + sum(P.*(alpha + g)./h_nu_xi)); % population of metastable level normalized by rho0
                               
                dP = u.*(alpha + g).*n2.*P... % medium gain
                    -u.*(alpha + edf.excess_loss*log(10)/10).*P... % attenuation
                    +ASEselect.*u.*g.*n2.*edf.Nmode*edf.h*edf.c./lamb*BWref; % ASE 
            end
            
            function res = bcfun(P0, PL, Pump, Signal, ASEf, ASEb)
                %% Calculate residual at boundaries
                
                % Boundary conditions for pump
                res = zeros(size(P0));
                if strcmpi(Pump.direction, 'forward')
                    res(1:Pump.N) = P0(1:Pump.N) - Pump.P.'; % P(z = 0) = Ppump
                else
                    res(1:Pump.N) = PL(1:Pump.N) - Pump.P.'; % P(z = L) = Ppump
                end
                
                % Boundary conditions for signal (always forward)
                idx = Pump.N + (1:Signal.N); % index signal
                res(idx) = P0(idx) - Signal.P.'; % P(z = 0) = Psignal
                
                % Forward ASE
                idx = idx + Signal.N; % index forward ASE
                res(idx) = P0(idx) - ASEf.P.'; % For single-stage, P(z = 0) = 0
                
                % Backward ASE
                idx = idx + Signal.N; % index backward ASE
                res(idx) = PL(idx) - ASEb.P.'; % For single-stage, P(z = L) = 0
            end
        end
        
        function [n2, z] = metastable_level_population(self, sol, Signal, Pump, ASE, verbose)
            %% Normalized metastable level population along the fiber (n2/nt)
            lamb = [Pump.wavelength, Signal.wavelength, ASE.wavelength, ASE.wavelength].'; % wavelength
            g = self.gain_coeff(lamb); % Gain coefficient
            alpha = self.absorption_coeff(lamb); % Absorption coefficient
                       
            % Photon energy times saturation parameter
            h_nu_xi = self.Ephoton(lamb).*self.sat_param(lamb); % h*nu*xi
                  
            z = sol.x;
            n2 = zeros(size(z));
            n2_approx = n2;
            for k = 1:size(sol.y, 2)
                P = sol.y(:, k);
                n2(k) = sum(P.*alpha./h_nu_xi)./(1 + sum(P.*(alpha + g)./h_nu_xi)); % population of metastable level normalized by rho0
                n2_approx(k) = alpha(1)/(alpha(1) + g(1)); % considering full inversion
            end
            
            if exist('verbose', 'var') && verbose
                figure(103), box on, hold on
                plot(z, n2, 'DisplayName', 'Numerical')
                plot(z, n2_approx, 'DisplayName', 'Full inversion approximation')
                xlabel('Distance (m)')
                ylabel('Normalized metastable level population')
                legend('-DynamicLegend')
            end  
        end
                
        %% EDF properties
        function [Gamma, confinement_factor] = overlap(self, lamb)
            %% Mode/doping-region overlap integral
            confinement_factor = self.doping_radius./self.mode_radius(lamb);
            
            Gamma = 1 - exp(-(confinement_factor).^2); 
            % follows from [3, eq. (1.220)] assuming step-like doping and
            % Gaussian envelope approximation for the mode
        end
        
        function sat_param = sat_param(self, lamb)
            %% Saturation parameter as defined in [1, below eq. 20]
            % Also equal to pi*(doping_radius)^2*(Erbium density)/tau
            sat_param = self.Psat(lamb).*(self.absorption_coeff(lamb) + self.gain_coeff(lamb))./self.Ephoton(lamb);
        end
                   
        function Psat = Psat(self, lamb)
            %% Saturation power as defined in [1, eq. XX]
            Aeff = self.effective_area(lamb);
            Psat = self.Ephoton(lamb).*Aeff./(self.tau*(self.absorption_cross_sec(lamb) + self.emission_cross_sec(lamb)));
        end
        
        function E = Ephoton(self, lamb)
            %% Photon energy
            E = self.h*self.c./lamb;
        end
              
        function [wBessel, wGauss] = mode_radius(self, lamb)
            %% Compute mode radius using Bessel solution [2, Appendix B]
            % This approximation is called Gaussian envelope in [2].
            % The mode shape is still approximated as a Gaussian, but the
            % mode radius is equal to the actual 1/e radius obtained from
            % the Bessel solution.
            V = self.V(lamb);
            U = (1 + sqrt(2))*V./(1 + (4 + V.^4).^(0.25));
            W = sqrt(V.^2 - U.^2);
            wBessel = self.core_radius*((V.*besselk(1, W))./(U.*besselk(0, W))).*besselj(0, U);
            
            wGauss = (self.core_radius)*(0.65 + 0.1619*V.^(-1.5) + 2.879*V.^(-6));  % Agrawal, eq. 2.2.43
            % The Gaussian approximation is accurate to within 1% only for 1.2 < V < 2.4.
            % The Bessel solution is preferred [2, Appendix B]
        end
        
        function V = V(self, lamb)
            %% V parameter (normalized frequency)
            V = 2*pi./lamb*self.core_radius*self.NA; % Agrawal, eq. 2.2.33
        end
        
        function Aeff = effective_area(self, lamb)
            %% Effective area
            Aeff = pi*(self.mode_radius(lamb)).^2; % definition
        end
        
        function [alpha, alphadB] = absorption_coeff(self, lamb) 
            %% Absorption coefficient in 1/m and in dB/m
            if isfield(self.param, 'absorption_coeff_fun')
                alphadB = self.param.absorption_coeff_fun(lamb*1e9); % converts lamb to nm before calling function
                alphadB(lamb >= 970e-9 & lamb <= 990e-9) = 10/log(10)*self.cross_sec2coeff(self.abs_cs980nm, 980e-9); % assign cross-section for 980 nm directly, since absorption_coeff_fun is obtained around 1550nm
                if size(alphadB, 1) ~= size(lamb, 1)
                    alphadB = alphadB.'; % ensures that dimensions are consistent
                end
                alpha = alphadB*log(10)/10; % convert from dB/m to 1/m
            else % if absorption coefficient is not given in param, it must be calculated from absorption cross section
                alpha = self.cross_sec2coeff(self.absorption_cross_sec(lamb), lamb); 
                alphadB = 10*alpha/log(10); % convert from 1/m to dB/m
            end            
        end
        
        function [g, gdB] = gain_coeff(self, lamb) 
            %% Stimulated gain coefficient in 1/m and in dB/m
            if isfield(self.param, 'gain_coeff_fun')
                gdB = self.param.gain_coeff_fun(lamb*1e9); % converts lamb to nm before calling function
                gdB(lamb >= 970e-9 & lamb <= 990e-9) = 0; % gain coefficient near 980nm is assumed zero, so that two-level system can be used
                if size(gdB, 1) ~= size(lamb, 1)
                    gdB = gdB.'; % ensures that dimensions are consistent
                end
                g = gdB*log(10)/10; % convert from dB/m to 1/m
            else % if gain coefficient is not given in param, it must be calculated from emission cross section
                g = self.cross_sec2coeff(self.emission_cross_sec(lamb), lamb); % 1/m
                gdB = 10*g/log(10); % convert from 1/m to dB/m
            end 
        end
             
        function acs = absorption_cross_sec(self, lamb)
            %% Absorption cross section (m^2) evaluated at wavelength lamb
            if isfield(self.param, 'abs_cross_sec')
                acs = self.Gaussian_fit(lamb, self.param.abs_cross_sec);
                acs(lamb >= 970e-9 & lamb <= 990e-9) = self.abs_cs980nm; % assign cross-section for 980 nm directly, since cross-section curves in edf_selection.m are valid near 1550 nm only
            else % if abs_cross_sec not in parameters, absorption_cross_sec is not calculated
                acs = self.coeff2cross_sec(self.absorption_coeff(lamb), lamb);
            end
        end
        
        function ecs = emission_cross_sec(self, lamb)
            %% Emission cross section (m^2) evaluated at wavelength lamb
            if isfield(self.param, 'ems_cross_sec')
                ecs = self.Gaussian_fit(lamb, self.param.ems_cross_sec);
                ecs(lamb >= 970e-9 & lamb <= 990e-9) = 0; % emission cross-section near 980nm is assumed zero, so that two-level system can be used
            else % if ems_cross_sec not in parameters, calculate cross section from gain coefficient
                ecs = self.coeff2cross_sec(self.gain_coeff(lamb), lamb);
            end
        end
        
        function eta = cross_section_ratio(self, lamb)
            %% Ratio between emission and absorption cross sections
            eta = self.emission_cross_sec(lamb)./self.absorption_cross_sec(lamb); % eq 1.27 of Principles
        end
        
        function coeff = cross_sec2coeff(self, cs, lamb)
            %% Convert (emission/absorption) cross-section to (gain/absorption) coefficient
            rho0m = self.rho0*1e6; % Er concentration in m^3
            coeff = cs.*rho0m.*self.overlap(lamb); % [3, Eq. 1.224]
        end
        
        function cs = coeff2cross_sec(self, coeff, lamb)
            %% Convert (gain/absorption) coefficient to (emission/absorption) cross-section
            rho0m = self.rho0*1e6; % Er concentration in m^3
            cs = coeff./(rho0m.*self.overlap(lamb)); % [3, Eq. 1.224]
        end
            
        % Auxiliary methods
        function g = Gaussian_fit(~, lamb, p)
            %% Sum of Gaussians. The kth Gaussian is centered around l(k), with peak a(k), and FWHM of d(k)
            l = p.l;
            a = p.a;
            d = p.d;
            g = zeros(size(lamb));
            sig = d/2*sqrt(2*log(2)); % relation between variance and FWHM
            for k = 1:length(a)
                g = g + a(k)*exp(-(lamb-l(k)).^2/sig(k)^2);
            end
        end
        
        function plot(self, property, lamb)
            %% Plot |property| for wavelengths given in lamb
            % Input: 
            % - property: what to plot = 'cross_sections', 'coefficients'
            % - lamb (optional, default = 1.4um to 1.6um): wavelength
            if not(exist('lamb', 'var'))
                lamb = 1e-6*linspace(1.47, 1.58);
            end
            
            % Plot
            if strcmpi(property, 'cross-sections') || strcmpi(property, 'all')
                figure(101), hold on, box on
                plot(lamb*1e9, absorption_cross_sec(self, lamb))
                plot(lamb*1e9, emission_cross_sec(self, lamb))
                ylabel('Cross-section (m^2)')
                xlabel('Wavelength (nm)')
                legend('Absorption', 'Emission')
            end
            if strcmpi(property, 'coefficients') || strcmpi(property, 'all')     
                figure(102), hold on, box on
                [alpha, alphadB] = absorption_coeff(self, lamb);
                [g, gdB] = gain_coeff(self, lamb);
                plot(lamb*1e9, alpha)
                plot(lamb*1e9, g)
                ylabel('Coefficient (1/m)')
                legend('Absorption', 'Gain')
                xlabel('Wavelength (nm)')
%                 axis([1470 1570 0 5])
            end
            if strcmpi(property, 'saturation param') || strcmpi(property, 'all')     
                figure(103), box on
                plot(lamb*1e9, self.sat_param(lamb))
                ylabel('Saturation parameter')
                xlabel('Wavelength (nm)')
            end
            
            if strcmpi(property, 'mode radius')
                [wBessel, wGauss] = mode_radius(self, lamb);
                figure(104), box on, hold on
                plot(lamb*1e9, wBessel*1e6, lamb*1e9, wGauss*1e6)
                legend('Bessel solution', 'Gaussian approximation')
                xlabel('Wavelength (nm)')
                ylabel('Mode radius (\mu m)')
            end
        end
    end
end