classdef EDF
    %% Single-mode Erbium-doped fiber
    % EDF parameters are loaded from file using function load_parameters.m
      
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
        type='corning high na' % which fiber to use. 
        % Fibers available:  
        % - 'corning high na' (default): high NA EDF provided by Corning
        % - 'corning_type1': first fiber data  
        % - 'corning_exp' : fiber used in gain experiments
        % - {'giles_ge:silicate', 'giles_al:ge:silicate'}. Extracted from [1, Fig. 2].
        % - {'principles_type1', 'principles_type2', 'principles_type3'}. 3 types of Er3+ in alumino-germanosilicate glass. Figs. 4.20, 4.21, and 4.22 of [2]
        L % fiber length (m)
        excess_loss = 0 % (dB/m) excess loss due to splices for instance. 
        % Default values assume 'corning_edf' fiber
        gp_980nm = 0; % gain coefficient is assumed 0, so that two-level system can be used
        alphap_980nm = 4.272 % absorption cross section near 980 nm (dB/m). 
        core_radius = 1.64e-6 % Fiber core radius. e.g, 1.2 um in [1, Table 1], 1.4 um in [4, pg. 156]
        doping_radius = 1.38e-6 % Er3+ core radius. e.g., 1.2 um in [1, Table 1], 1.05um in [4, pg 156]
        rho0 = 5.51e18; % Er3+ concentraction (cm^3), e.g., 0.7e19 in [4, pg 156]
        NA = 0.23 % numerical aperture, e.g., 0.28 in [4, pg. 156]
        tau = 10e-3; % metastable lifetime in s        
        Nmode = 2 % number of modes (default = 2 for two polarizations)
        correlation_fit_file = 'corr_fun_fit.mat' % file containg fit of correlation function Gamma. Only used when modeling spectral hole burning
    end
         
    properties
        param % other physical parameters from fiber (loaded from load_parameters.m)
    end
        
    properties (Constant, Hidden)
        
        % Value obtainend from [4, pg 154]. Cross section curves in edf_select.m are only valid near 1550 nm
        maxL = 10; % maximum fiber length. Used to limit simulation
        h = 6.62606957e-34; % Planck
        q = 1.60217657e-19; % electron charge
        c = 299792458;      % speed of light
    end
    
    methods
        function obj = EDF(L, type)
            %% Constructor
            %       
            obj.L = L;
            if exist('type', 'var')
                obj.type = type;
            end
            obj = load_parameters(obj);
        end
        
        % set methos
        function self = set.type(self, t)
            self.type = t;
            self = load_parameters(self);
        end
                
        %% Semi-analytical model
        function [GaindB, Psignal_out, Ppump_out, dGaindB, dGaindB_L] = semi_analytical_gain(self, Pump, Signal)
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
            % Outputs:
            % - GaindB: gain of each Signal channel in dB
            % - Psignal_out: output power of each signal in W
            % - Ppump_out: 
            % - dGaindB: gradient of the gain in dB with respect to the signal power in W.
            % dGaindB is a N x N matrix where dGain(i, j) = partial Gain_i /
            % partial P_j
            % - dGaindB_L: gradient of the gain in dB with respect to the
            % EDF length
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
            a = (alpha + g)./xi;
            b = alpha*self.L;
            [Qout, ~, exitflag] = fzero(@(Qout) Qout - sum(Qin_k.*exp(a*(Qin - Qout) - b)), Qout0);
            
            if exitflag ~= 1
                warning('EDF/semi_analytical_gain: simulation exited with exitflag = %d\n', exitflag)
            end
            
            % Calculate the output flux for each individual signal & pump
            % using equation (24) of [1]
            Gain = exp(a*(Qin - Qout) - b);
            Qout_k = Qin_k.*Gain;
            
            Ppump_out = Qout_k(1:Pump.N).*self.Ephoton(Pump.wavelength);
            Psignal_out = Qout_k(Pump.N+1:end).*self.Ephoton(Signal.wavelength);
            
            GaindB = 10*log10(Gain(Pump.N+1:end));
            GaindB = GaindB - self.L*self.excess_loss; % acounts for excess loss
            
            if nargout >= 4 % gradient was requested
                dGain = (a(Pump.N+1:end))/(1 + sum(Qout_k.*a)); % 1 x N vector 
                dGain = dGain.*((1 - Gain(Pump.N+1:end))./self.Ephoton(Signal.wavelength)).'; % N X N matrix dGain(i,j) = partial Gain_i / partial P_j 
                dGaindB = 10/log(10)*dGain; % gain cancelled out
                
                dQout_L = -sum(Qout_k.*alpha)/(1 + sum(Qout_k.*a));
                dGaindB_L = -10/log(10)*(alpha(Pump.N+1:end) + a(Pump.N+1:end)*dQout_L).'; 
            end
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
            %% Analytical excess noise [1, eq. (31)]
            % Assumes that the EDF is fully inverted
            gp = self.gain_coeff(Pump.wavelength);
            alphap = self.absorption_coeff(Pump.wavelength);
            
            gs = self.gain_coeff(Signal.wavelength);
            alphas = self.absorption_coeff(Signal.wavelength);
            
            nsp = 1./(1 - gp*alphas./(alphap*gs));
        end
        
        function NFdB = analytical_noise_figure(self, Pump, Signal)
            %% Analytical noise figure in dB
            nsp = analytical_excess_noise(self, Pump, Signal);
            GdB = semi_analytical_gain(self, Pump, Signal);
            G = 10.^(GdB/10);
            NFdB = 10*log10(2*nsp.*(G-1)./G); % By definition. Note that for high gain, NF is approximately 2nsp
        end   

        function Pase = analytical_ASE_PSD(self, Pump, Signal, nsp_correction)
            %% Analytical expression for ASE PSD [1, eq. (32)]
            % The optional parameter nsp_correction is a multiplicative
            % factor to correct the analytical nsp. The analytical nsp
            % assumes that the EDF is uniformly inverted, which leads to
            % more optimistic values of nsp
            % Usually nsp_correction = 1.2 for 980nm and nsp_correction =
            % 1.5 for 1480 nm leads to good results
            if not(exist('nsp_correction', 'var'))
                nsp_correction = 1;
            end
            
            nsp = nsp_correction*analytical_excess_noise(self, Pump, Signal);
            GdB = semi_analytical_gain(self, Pump, Signal);
            Pase = 2*nsp.*(10.^(GdB/10)-1).*Signal.Ephoton;
            Pase = max(Pase, 0); % non-negativity constraint
        end
               
        function Lopt = optimal_length(self, Pump, Signal, spanAttdB)
            %% Optimize EDF length to maximize number of channels with non-zero power that have gain above spanAttdB
            % This function uses the semi_analytical_gain() to obtain the gain
            if length(spanAttdB) == 1
                spanAttdB = spanAttdB*ones(size(Signal.P));
            end
            Temp = self;
            [Lopt, ~, exitflag] = fminbnd(@(L) objective(Temp, L, spanAttdB), 0, self.maxL);
           
            if exitflag ~= 1
                warning('EDFA/optimal_length: could not solve for Qout. Simulation exited with exitflag = %d\n', exitflag)
            end
            
            function y = objective(Temp, L, spanAttdB)
                Temp.L = L;
                GaindB = Temp.semi_analytical_gain(Pump, Signal); 
                y = -sum(GaindB(Signal.P ~= 0) >= spanAttdB(Signal.P ~= 0)); 
            end
        end
        
        %% Numerical models 
        function Gamma = calc_corr_fun(self, lamb)
            %% Compute correlation function at wavelenghts lamb
            lamb = lamb*1e9;
            [X, Y] = meshgrid(lamb, lamb);

            CorrFun = load(self.correlation_fit_file);
            Gamma = -5*CorrFun.Gamma_fit(X, Y);
%             Gamma = Gamma*x(1) + x(2);
            Gamma(:, 1) = 0;
            Gamma(1, :) = 0; % erase pump correlation coefficients i.e., pump does not induce SHB of the signal
        end
        
        function [GaindB, Ppump_out, Psignal_out, Pase, sol]...
                = propagate(self, Pump, Signal, ASEf, ASEb, BWref, model, Nmesh, verbose)
            %% Calculate Gain and noise PSD by solving coupled nonlinear first-order differential equations that follow from the SCD model                                             
            % Inputs:
            % - Pump: instance of class Channels corresponding to Pump
            % - Signal: instance of class Channels corresponding to Signal
            % - ASEf: instance of class Channels corresponding to forward ASE
            % - ASEb: instance of class Channels corresponding to backward ASE
            % - BWref: reference bandwidth for computing ASE power
            % - system (optional, default='two-level'): laser system to simulate either 'two-level' or 'three-level'
            % - Nmesh (optional, default=50e3/length(lamb)): number of mesh points to use when solving BVP
            % - verbose (optiona, default=false): whether to plot results
            
            lamb = [Pump.wavelength, Signal.wavelength, ASEf.wavelength, ASEb.wavelength].'; % wavelength
            u = [Pump.u, Signal.u, ASEf.u, ASEb.u].'; % propation direction
            g = self.gain_coeff(lamb); % Gain coefficient
            alpha = self.absorption_coeff(lamb); % Absorption coefficient
            
            % ASE term is only included in the forward and backward ASE channels
            ASEselect = ones(3*Signal.N, 1);
            ASEselect(1:Signal.N) = 0;
            ASEselect = logical(ASEselect);
            
            % Photon energy times saturation parameter
            xi = self.sat_param(lamb);
            h_nu_xi = self.Ephoton(lamb).*xi; % h*nu*xi
            
            % Solver
            if not(exist('Nmesh', 'var'))
                Nmesh = floor(50e3/length(lamb));
                % Maximum number of mesh points allowed when solving the BVP
            end
                
            z = [0, self.L];
            P0 = [Pump.P Signal.P ASEf.P ASEb.P].';
            solinit = bvpinit(z, P0); % initial guess
            
            if not(exist('model', 'var'))
                model = 'two-level';
            end
            
            switch lower(model)
                case 'three-level' % three-level system
                    options = bvpset('Vectorized', 'off', 'NMax', Nmesh);
                    sol = bvp4c(@(z, P) three_level_system(z, P, self, lamb, h_nu_xi, u, g, alpha, BWref),... % differential equation
                        @(P0, PL) bcfun(P0, PL, Pump, Signal, ASEf, ASEb),... % boundary conditions
                        solinit, options); % initial guess
                case 'two-level' % two-level system
                    options = bvpset('Vectorized', 'off', 'NMax', Nmesh);
                    sol = bvp4c(@(z, P) two_level_system(z, P, self, lamb, h_nu_xi, u, g, alpha, BWref),... % differential equation
                        @(P0, PL) bcfun(P0, PL, Pump, Signal, ASEf, ASEb),... % boundary conditions
                        solinit, options); % initial guess
                case 'three-level+shb' % three-level system including spectral hole burning
                    Gamma = self.calc_corr_fun(lamb);
                    options = bvpset('Vectorized', 'off', 'NMax', Nmesh);
                    sol = bvp4c(@(z, P) three_level_system_SHB(z, P, self, lamb, h_nu_xi, u, g, alpha, BWref),... % differential equation
                        @(P0, PL) bcfun(P0, PL, Pump, Signal, ASEf, ASEb),... % boundary conditions
                        solinit, options); % initial guess
                case 'two-level+shb' % two-level system including spectral hole burning
                    Gamma = self.calc_corr_fun(lamb);
                    options = bvpset('Vectorized', 'off', 'NMax', Nmesh);
                    sol = bvp4c(@(z, P) two_level_system_SHB(z, P, self, lamb, h_nu_xi, u, g, alpha, BWref),... % differential equation
                        @(P0, PL) bcfun(P0, PL, Pump, Signal, ASEf, ASEb),... % boundary conditions
                        solinit, options); % initial guess
                otherwise
                    error('edf.propagate: unknown model type')
            end
                
            if any(sol.y(:) < 0)
                warning('EDF/propagate: solution contains negative power')
            end
 
            if strcmpi(Pump.direction, 'forward')
                Ppump_out = sol.y(1:Pump.N, end).';
            else
                Ppump_out = sol.y(1:Pump.N, 1).';
            end
                
            Psignal_out = sol.y(Pump.N + (1:Signal.N), end).';
            Pase = sol.y(Pump.N + Signal.N + (1:Signal.N), end).'; % only forward ASE
            Pase_backward = sol.y(Pump.N + 2*Signal.N + (1:Signal.N), 1).';
            GaindB = 10*log10(Psignal_out./Signal.P); 
            
            if exist('verbose', 'var') && verbose
                figure(243)
                subplot(221), hold on, box on % pump evolution
                plot(sol.x, 1e3*sol.y(1:Pump.N, :))
                xlabel('z (m)')
                ylabel('Pump power (mW)')
                xlim([0 self.L])
                title('Pump evolution')
                
                subplot(222), hold on, box on % signal evolution
                plot(sol.x, Watt2dBm(sol.y(Pump.N + (1:Signal.N), :)) - Signal.PdBm.') % broadcast
                xlabel('z (m)')
                ylabel('Signal gain (dB)')
                xlim([0 self.L])
                title('Signal gain evolution')
                
                subplot(223), hold on, box on % ASEf evolution
                plot(sol.x, Watt2dBm(sol.y(Pump.N + Signal.N + (1:Signal.N), :)))
                xlabel('z (m)')
                ylabel('ASE power (dBm)')
                xlim([0 self.L])
                title('Forward ASE evolution')
                
                subplot(224), hold on, box on % ASEb evolution
                plot(sol.x, Watt2dBm(sol.y(Pump.N + 2*Signal.N + (1:Signal.N), :)))
                xlabel('z (m)')
                ylabel('ASE power (dBm)')
                xlim([0 self.L])
                title('Backward ASE evolution')
                
                %
                figure(244)
                subplot(311), hold on, box on % gain
                plot(Signal.wavelength*1e9, GaindB)
                xlabel('Wavelength (nm)')
                ylabel('Gain (dB)')
                title('Gain')
                
                subplot(312), hold on, box on % Pout
                PoutdBm = Signal.PdBm + GaindB;
                plot(Signal.wavelength*1e9, PoutdBm)
                xlabel('Wavelength (nm)')
                ylabel('Output signal power (dB)')
                title('Signal output')
                ylim([floor(min(PoutdBm(PoutdBm > -50))) ceil(max(PoutdBm))])
                
                subplot(313), hold on, box on % ASE
                plot(Signal.wavelength*1e9, Watt2dBm(Pase), 'DisplayName', 'Forward')
                plot(Signal.wavelength*1e9, Watt2dBm(Pase_backward), 'DisplayName', 'Backward')
                legend('-dynamiclegend', 'Location', 'Best')
                xlabel('Wavelength (nm)')
                ylabel('ASE (dBm)')
                title('ASE') 
            end
            
            function dP = three_level_system(~, P, edf, lamb, h_nu_xi, u, g, alpha, BWref)
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
                p = 1; % pump
                s = 2:length(lamb); % signal & ASE
                Ephoton = edf.h*edf.c./lamb(s);

                ap = alpha(p)*(1 + sum(P(s).*g(s)./h_nu_xi(s)))./(1 + sum(P.*(alpha + g)./h_nu_xi));
                ep = 0;
                es = g(s).*(P(p)*(alpha(p) + g(p))/h_nu_xi(p) + sum(alpha(s).*P(s)./h_nu_xi(s)))./(1 + sum(P.*(alpha + g)./h_nu_xi));
                as = alpha(s).*(1 + sum(g(s).*P(s)./h_nu_xi(s)))./(1 + sum(P.*(alpha + g)./h_nu_xi));

                dPp = u(p)*(ep - ap - edf.excess_loss*log(10)/10)*P(p); % Pump
                dPs = u(s).*(es - as - edf.excess_loss*log(10)/10).*P(s); % Signal
                dPs(ASEselect) = dPs(ASEselect) + 2*u(s(ASEselect)).*es(ASEselect).*(BWref*Ephoton(ASEselect)); % Include ASE term in ASE channels

                dP = [dPp; dPs];
            end

            function dP = two_level_system(~, P, edf, lamb, h_nu_xi, u, g, alpha, BWref)
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
                    +double([zeros(Pump.N, 1); ASEselect]).*u.*g.*n2.*edf.Nmode*edf.h*edf.c./lamb*BWref; % ASE 
                
                1;
            end
           
            function dP = three_level_system_SHB(~, P, edf, lamb, h_nu_xi, u, g, alpha, BWref)
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
                               
                Rup = sum(alpha.*P);
                Rdown = sum(g.*P + xi.*edf.h*edf.c./lamb);
                Rup_prime = Gamma*(alpha.*P);
                Rdown_prime = Gamma*(g.*P);
                
                at = alpha.*(1 + (Rup*Rdown_prime - Rdown*Rup_prime)./(Rdown*(Rup+Rdown)));
                gt = g.*(1 - (Rup*Rdown_prime - Rdown*Rup_prime)./(Rup*(Rup+Rdown)));
                
                %
                p = 1; % pump
                s = 2:length(lamb); % signal & ASE
                Ephoton = edf.h*edf.c./lamb(s);
                                
                ap = at(p)*(1 + sum(P(s).*gt(s)./h_nu_xi(s)))./(1 + sum(P.*(at + gt)./h_nu_xi));
                ep = 0;
                es = gt(s).*(P(p)*(at(p) + gt(p))/h_nu_xi(p) + sum(at(s).*P(s)./h_nu_xi(s)))./(1 + sum(P.*(at + gt)./h_nu_xi));
                as = at(s).*(1 + sum(gt(s).*P(s)./h_nu_xi(s)))./(1 + sum(P.*(at + gt)./h_nu_xi));
                                    
                dPp = u(p)*(ep - ap - edf.excess_loss*log(10)/10)*P(p); % Pump
                dPs = u(s).*(es - as - edf.excess_loss*log(10)/10).*P(s); % Signal
                dPs(ASEselect) = dPs(ASEselect) + 2*u(s(ASEselect)).*es(ASEselect).*(BWref*Ephoton(ASEselect)); % Include ASE term in ASE channels
                
                dP = [dPp; dPs];
           end
            
           function dP = two_level_system_SHB(~, P, edf, lamb, h_nu_xi, u, g, alpha, BWref)
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
                               
                Rup = sum(alpha.*P);
                Rdown = sum(g.*P + xi.*edf.h*edf.c./lamb);
                Rup_prime = Gamma*(alpha.*P);
                Rdown_prime = Gamma*(g.*P);
                
                at = alpha.*(1 + (Rup*Rdown_prime - Rdown*Rup_prime)./(Rdown*(Rup+Rdown)));
                gt = g.*(1 - (Rup*Rdown_prime - Rdown*Rup_prime)./(Rup*(Rup+Rdown)));
                
                n2 = sum(P.*at./h_nu_xi)./(1 + sum(P.*(at + gt)./h_nu_xi)); % population of metastable level normalized by rho0
                    
                dP = u.*(at + gt).*n2.*P... % medium gain
                    -u.*(at + edf.excess_loss*log(10)/10).*P... % attenuation
                    +double([0; ASEselect]).*u.*gt.*n2.*edf.Nmode*edf.h*edf.c./lamb*BWref; % ASE 
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
           
        function [GaindB, Ppump_out, Psignal_out, Pase, sol]...
                = two_level_system(varargin)
            %[GaindB, Ppump_out, Psignal_out, Pase, sol] = two_level_system(self, Pump, Signal, ASEf, ASEb, BWref, Nmesh, verbose)
            %% Legacy; use method "propagate" instead
            
            if nargin > 6
                v = [varargin(1:6), 'two-level', varargin(7:nargin)];
            else
                v = [varargin(1:6) 'two-level'];
            end
            1;
            [GaindB, Ppump_out, Psignal_out, Pase, sol] = propagate(v{:});
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
                figure(203), box on, hold on
                plot(z, n2, 'DisplayName', 'Numerical')
                plot(z, n2_approx, 'DisplayName', 'Full inversion approximation')
                xlabel('Distance (m)')
                ylabel('Normalized metastable level population')
                legend('-DynamicLegend')
            end  
        end
        
        function [nsp, NFdB] = excess_noise(self, Pump, Signal, verbose)
            %% Compute excess noise factor (nsp) numerically [Principles, eq 5.31]
            % EDFA gain and forward ASE should be measured in a narrow
            % linewidth (e.g., 1nm)
            BWref = 50e9; 
            
            Signal.P(Signal.P == 0) = 1e-8; % set channels with small power
            ASEf = Channels(Signal.wavelength, 0, 'forward');
            ASEb = Channels(Signal.wavelength, 0, 'backward');
%             
            [GaindB, ~, ~, Pase] = self.propagate(Pump, Signal, ASEf, ASEb, BWref, 'three-level', 100, verbose);
                        
            % Compute excess noise
            G = 10.^(GaindB/10);            
            nsp = 1./(G-1).*Pase./(2*BWref*self.h*self.c./Signal.wavelength);
            NFdB = 10*log10(2*nsp.*(G-1)./G); % [Principles, eq. 5.34]
            
            if exist('verbose', 'var') && verbose
                nsp_analytical = self.analytical_excess_noise(Pump, Signal);
                NFdB_analytical = self.analytical_noise_figure(Pump, Signal);
                figure(187)
                subplot(211), hold on, box on
                plot(Signal.lnm, nsp, 'DisplayName', 'Numerical')
                plot(Signal.lnm, nsp_analytical, 'DisplayName', 'Anlytical')
                xlabel('Wavelength (nm)')
                ylabel('Excess noise factor')
                legend('-dynamiclegend')
                
                subplot(212), hold on, box on
                plot(Signal.lnm, NFdB, 'DisplayName', 'Numerical')
                plot(Signal.lnm, NFdB_analytical, 'DisplayName', 'Anlytical')
                xlabel('Wavelength (nm)')
                ylabel('Noise Figure (dB)')
                legend('-dynamiclegend')
            end
        end
        
        function [PCE, PCEmax] = power_conversion_efficiency(~, Pump, SignalIn, SignalOut)
            %% Power conversion efficiency            
            PCE = sum(SignalOut.P - SignalIn.P)/Pump.P; % definition [Principles, eq. 5.22]
            PCEmax = Pump.wavelength/min(SignalIn.wavelength(SignalIn.P ~= 0)); % theoretical maximum by conservation of energy [Principles, 5.23]
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
            sat_param = pi*(self.doping_radius)^2*(1e6*self.rho0)/self.tau;
%             sat_param = 3.5e15; % experimental value for Corning fiber
%             "Corning (NEW)"
            % Below is another form of calculating sat_param when ignoring
            % the transverse variation of the mode intensity across the
            % doping profile
%             sat_param = self.Psat(lamb).*(self.absorption_coeff(lamb) + self.gain_coeff(lamb))./self.Ephoton(lamb);
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
            if any(lamb > 1570)
                warnig('Abosorption and gain coefficients can only be calculated confidently for wavelengths below 1570 nm.')
            end
            if isfield(self.param, 'absorption_coeff_fun')
                alphadB = self.param.absorption_coeff_fun(lamb*1e9); % converts lamb to nm before calling function
                if isfield(self.param, 'pump_absorption_coeff_fun')
                    alphadB(lamb <= 990e-9) = self.param.pump_absorption_coeff_fun(1e9*lamb(lamb <= 990e-9)); % converts lamb to nm before calling function
                else
                    alphadB(lamb <= 990e-9) = self.alphap_980nm; % assign cross-section for 980 nm directly, since absorption_coeff_fun is obtained around 1550nm
                end
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
            if any(lamb > 1570)
                warnig('Abosorption and gain coefficients can only be calculated confidently for wavelengths below 1570 nm.')
            end
            if isfield(self.param, 'gain_coeff_fun')
                gdB = self.param.gain_coeff_fun(lamb*1e9); % converts lamb to nm before calling function
                if isfield(self.param, 'pump_gain_coeff_fun')
                    gdB(lamb <= 990e-9) = self.param.pump_gain_coeff_fun(1e9*lamb(lamb <= 990e-9)); % converts lamb to nm before calling function
                else
                    gdB(lamb <= 990e-9) = self.gp_980nm; % assign cross-section for 980 nm directly
                end
                
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
                acs(lamb >= 970e-9 & lamb <= 990e-9) = self.coeff2cross_sec(log(10)/10*self.alphap_980nm, 980e-9); % assign cross-section for 980 nm directly, since cross-section curves in edf_selection.m are valid near 1550 nm only
            else % if abs_cross_sec not in parameters, absorption_cross_sec is not calculated
                acs = self.coeff2cross_sec(self.absorption_coeff(lamb), lamb);
            end
        end
        
        function ecs = emission_cross_sec(self, lamb)
            %% Emission cross section (m^2) evaluated at wavelength lamb
            if isfield(self.param, 'ems_cross_sec')
                ecs = self.Gaussian_fit(lamb, self.param.ems_cross_sec);
                ecs(lamb >= 970e-9 & lamb <= 990e-9) = self.coeff2cross_sec(log(10)/10*self.gp_980nm, 980e-9); % emission cross-section near 980nm is assumed zero, so that two-level system can be used
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
                g = g + a(k)*exp(-((lamb-l(k))/sig(k)).^2);
            end
        end
        
        function plot(self, property, lamb)
            %% Plot |property| for wavelengths given in lamb
            % Input: 
            % - property: what to plot = 'cross_sections', 'coefficients'
            % - lamb (optional, default = 1.47um to 1.57um): wavelength
            if not(exist('lamb', 'var'))
                lamb = 1e-6*linspace(1.47, 1.57);
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
            
            if strcmpi(property, 'mode radius')
                [wBessel, wGauss] = mode_radius(self, lamb);
                figure(103), box on, hold on
                plot(lamb*1e9, wBessel*1e6, lamb*1e9, wGauss*1e6)
                legend('Bessel solution', 'Gaussian approximation')
                xlabel('Wavelength (nm)')
                ylabel('Mode radius (\mu m)')
            end
        end
    end
end
