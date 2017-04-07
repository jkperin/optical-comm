%% Avalanche photodiode
classdef apd
    properties
        Gain % Gain
        ka   % impact ionization factor
        BW0  % Low-gain bandwidth
        GainBW  % Gain Bandwidth product
        R    % responsivity
        Id   % dark current
    end
    
    properties (Dependent)
        Fa % excess noise factor
        Geff % effective gain = Gain*Responsivity
        BW   % Bandwidth
        GaindB % Gain in dB 
    end
    
    properties (Constant, Hidden)
        h = 6.62606957e-34; % Planck
        q = 1.60217657e-19; % electron charge
        c = 299792458;      % speed of light
    end
    
    properties (Constant, GetAccess=private, Hidden)
        cdf_accuracy = 1-1e-4; % Required accuracy of the cdf (i.e., pmf will have the minimum number of points that satistifes that sum pmf > cdf_accuracy.)
        Niterations = 1e6; % maximum number of iterations in a while loop
        Ptail = 1e-6; % probability of clipped tail
    end
    
    properties (Dependent, GetAccess=private, Hidden)
        a % auxiliary variable G = 1/(1-ab)
        b % auxiliary variable b = 1/(1-ka)
    end
   
    methods
        function this = apd(GaindB, ka, BW, R, Id)
            %% Class constructor
            % Input:
            % - GaindB = gain in dB
            % - ka = impact ionization factor
            % - BW (optional, default = Inf) = if number, then BW specifies bandwidth
            %   if BW is 2x1 vector, then 1st element is the low-gain bandwidth 
            %   and the second is the gain-bandwidth product    
            % - R (optional, default = 1 A/W) = responsivity
            % - Id (optional defualt = 0 nA) = dark current (A)
            
            % Required parameters: GaindB and ka
            this.Gain = 10^(GaindB/10);
            this.ka = ka;
            
            if exist('BW', 'var')
                if length(BW) == 1
                    this.BW0 = BW;
                    this.GainBW = Inf;
                else
                    this.BW0 = BW(1);
                    this.GainBW = BW(2);
                end
            else
                this.BW0 = Inf;
                this.GainBW = Inf;
            end
                        
            if exist('R', 'var')
                this.R = R;
            else 
                this.R = 1;
            end
            
            if exist('Id', 'var')
                this.Id = Id;
            else 
                this.Id = 0;
            end
        end
        
        function APDtable = summary(self)
            %% Generate table summarizing class values
            disp('-- APD class parameters summary:')
            rows = {'Gain'; 'Impact ionization factor'; 'Low-gain bandwidth';...
                'Gain-bandwidth product'; 'Responsivity'; 'Dark current'};
            Variables = {'Gain'; 'ka'; 'BW0'; 'GainBW'; 'R'; 'Id'};
            Values = [self.Gain; self.ka; self.BW0/1e9; self.GainBW/1e9; self.R; self.Id*1e9];
            Units = {''; ''; 'GHz'; 'GHz'; 'A/W'; 'nA'};

            APDtable = table(Variables, Values, Units, 'RowNames', rows);
        end
        
        %% Main Methods
        function l = lambda(this, P, dt)
            %% Rate of Possion process for a given P in a interval dt
            % From Personick, "Statistics of a General Class of Avalanche Detectors With Applications to Optical Communication"
            l = (this.R*P + this.Id)*dt/this.q; 
        end
        
        function Hapd = H(this, f)
            %% APD frequency response. Normalized to have unit gain at DC
            if isinf(this.BW)
                Hapd = ones(size(f));
            else
%                 Hapd = 1./sqrt(1 + (f/this.BW).^2);
                Hapd = 1./(1 + 1j*(f/this.BW))...
                    .*exp(1j*2*pi*f.*(this.BW./(2*pi*(this.BW^2 + f.^2)))); % remove group delay
            end
        end
        
        function hapd = ht(this, t)
            %% APD impulse response
            if isinf(this.BW)
                hapd = double(t == 0);
            else
%                 Hapd = 1./sqrt(1 + (f/this.BW).^2);
                hapd = this.BW*exp(-t*this.BW);
                hapd(t < 0) = 0;
            end            
        end
        
        function [Hw, y] = Hwhitening(self, f, P, N0, x)
            %% Design whitening filter
            % Inputs:
            % - f = frequency vector
            % - P = power at the APD input
            % - N0 = power spectral density of thermal noise
            % - x (optional) if provided y = x filtered by Hw
            if ~isinf(self.BW) % Only included if APD bandwidth is not inf
                r = N0/self.varShot(P, 1); % Ratio between thermal and shot noise PSD
                % Note: in calculating shot noise PSD the average power is used
                Hw = sqrt((1 + r)./(r + abs(self.H(f)).^2));
                [Hw, groupdelay] = Hgrpdelay(Hw, f);

                if exist('x', 'var')
                    y = ifft(fft(x).*ifftshift(Hw));
                end
            else
                Hw = 1;
                y = [];
                if exist('x', 'var')
                    y = x;
                end
            end
        end
        
        function bw = noisebw(this) 
            %% Noise bandwidth (one-sided)
            bw = pi/2*this.BW; % pi/2 appears because this.Hapd energy is pi*BW
        end
        
        function sig2 = varShot(this, Pin, Df)
            %% Shot noise variance
            % Inputs:
            % - Pin = input power of photodetector (W)
            % - Df = noise bandwidth (Hz)
            sig2 = 2*this.q*this.Gain^2*this.Fa*(this.R*Pin + this.Id)*Df; % Agrawal 4.4.17 (4th edition)
        end

        function output = detect(this, Ein, fs, noise_stats, N0)
            %% Direct detection
            % Inputs:
            % - Ein = input electric field
            % - fs = sampling frequency (Hz)
            % - noise_stats = 'gaussian', 'doubly-stochastic' (not
            % implemented), or 'no noise'
            % - N0 (optional, if provided, thermal noise of psd N0 is added
            % after direct detection) = thermal noise psd      
            
            if any(size(Ein) == 2) % two pols
                if size(Ein, 1) ~= 2
                    Ein = Ein.'; % put in 2 x N format
                end
                Pin = sum(abs(Ein).^2, 1);
            else
                Pin = abs(Ein).^2;
            end     
            
            switch noise_stats 
                case 'gaussian'
                    % Assuming Gaussian statistics    
                    output = this.R*this.Gain*Pin + sqrt(this.varShot(Pin, fs/2)).*randn(size(Pin));
                    
                case 'doubly-stochastic' % DEPRECTED
                    error('apd/detect: noise_stats = doubly-stochastic is deprected.');
%                     % uses saddlepoint approximation to obtain pmf in order to generate 
%                     % output distributed according to that pmf
%                     Plevels = unique(Pin);
% 
%                     output = zeros(size(Pin));
%                     for k = 1:length(Plevels)
%                         [px, x] = this.output_pdf_saddlepoint(Plevels(k), fs, 0); % doesn't include thermal noise here
% 
%                         cdf = cumtrapz(x, px);
% 
%                         pos = (Pin == Plevels(k));
% 
%                         % Sample according to pmf px
%                         u = rand(sum(pos), 1); % uniformly-distributed
%                         dist = abs(bsxfun(@minus, u, cdf));
%                         [~, ix] = min(dist, [], 2);
%                         output(pos) = x(ix); % distributed accordingly to px
%                     end
                case 'no noise'
                    % Only amplifies and filters the signal (no noise).
                    % This is used in estimating the BER and debugging
                    output = this.R*this.Gain*Pin;
                    
                otherwise 
                    error('apd/detect: Invalid Option!')
            end
            
            % Frequency
            df = fs/length(Pin);
            f = (-fs/2:df:fs/2-df);
            if size(output, 1) > size(output, 2)
                f = f.';
            end
            
            % APD frequency response
            if ~isinf(this.BW)
                output = real(ifft(fft(output).*ifftshift(this.H(f))));
                % H has unit gain at DC
            end
                        
            % Add thermal noise if N0 was provided
            if exist('N0', 'var')
                output = output + sqrt(N0*fs/2).*randn(size(Pin)); % includes thermal noise
            end 
        end
              
        function noise_std = stdNoise(this, Hrx, Hff, N0, RIN, sim)
            %% Returns function handle noise_std, which calculates the noise standard deviation for a given power level P. 
            % !! The power level P is assumed to be referred to after the APD.
            % Thermal, shot and RIN are assumed to be white AWGN.
            % Inputs:
            % - Hrx = receiver filter evaluated at sim.f e.g., whitening filter and matched filter 
            % - Hff = equalizer frequency response evaluated at sim.f
            % - N0 = thermal noise PSD
            % - RIN = RIN in dB/Hz (if empty RIN is not included)
            % - sim = sim struct            
            
            Htot = Hrx.*Hff; 
            Df  = 1/2*trapz(sim.f, abs(Htot).^2); % filter noise BW (includes noise enhancement penalty)
            if ~isinf(this.BW)
                Dfshot = 1/2*trapz(sim.f, abs(this.H(sim.f).*Htot).^2);
            else
                Dfshot = Df;
            end
            DfRIN = Dfshot; % !! approximated: needs to include fiber response
            
            %% Noise calculations
            % Thermal noise
            varTherm = N0*Df; % variance of thermal noise

            % RIN
            if ~isempty(RIN) && isfield(sim, 'RIN') && sim.RIN
                varRIN =  @(Plevel) 10^(RIN/10)*Plevel.^2*DfRIN;
            else
                varRIN = @(Plevel) 0;
            end

            % Shot
            varShot = @(Plevel) this.varShot(Plevel/(this.Gain*this.R), Dfshot);
            % Note: Plevel is divided by APD gain to obtain power at the
            % apd input, which determines the noise variance
            
            % Noise std for the level Plevel
            noise_std = @(Plevel) sqrt(varTherm + varRIN(Plevel) + varShot(Plevel));
        end
        
        function [Gopt, mpam] = optGain(this, mpam, tx, fiber, rx, sim)
            %% Optimize APD gain: Given target BER finds APD gain that leads to minimum required optical power
            disp('Optimizing APD gain for sensitivity');
                   
            % Find Gapd that minimizes required power to achieve target BER
            if mpam.optimize_level_spacing
                % Level spacing optimization ensures that target BER is met for a given gain
                % Thus, find APD gain that leads to minimum average
                % PAM levels
                [Gopt, ~, exitflag] = fminbnd(@(Gapd) ...
                    this.optimize_PAM_levels(Gapd, mpam, tx, fiber, rx, sim), eps, maxGain(this, mpam.Rs/5));  

                [~, mpam] = this.optimize_PAM_levels(Gopt, mpam, tx, fiber, rx, sim);
            else
                % Adjust power to get to target BER
                [Gopt, ~, exitflag] = fminbnd(@(Gapd) fzero(@(PtxdBm)...
                    log10(this.calc_apd_ber(PtxdBm, Gapd, mpam, tx, fiber, rx, sim)) - log10(sim.BERtarget), -20), 1, maxGain(this, mpam.Rs/5));
            end

            % Check whether solution is valid
            if exitflag ~= 1
                warning('apd/optGain: APD gain optimization did not converge (exitflag = %d)\n', exitflag);
            end 

            assert(Gopt >= 0, 'apd/optGain: Negative gain found while optimizing APD gain')
                    
            % Auxiliary function
            function Gmax = maxGain(apd, minBW)
                % Max gain allowed during gain optimization
                if isinf(apd.GainBW)
                    Gmax = 100;
                else
                    Gmax = apd.GainBW/minBW; 
                end
            end
        end 
        
        function [Pmean, mpam] = optimize_PAM_levels(this, Gapd, mpam, Tx, Fiber, Rx, sim)
            %% Calculate optimal level spacing for a given APD gain
            % Noise whitening filter depends on optical power, so the levels
            % must be calculated iteratively.
            this.Gain = Gapd;
                       
            %% Iterate until convergence
            Pmean = [0 1];
            tol = 1e-6;
            n = 1;
            maxIterations = 50;
            Tx.Ptx = 1e-3; % initial power 
            Fiber.att = @(l) 0; % disregard attenuation
            Pdiff = Inf;
            mpam = mpam.adjust_levels(Tx.Ptx, Tx.Mod.rexdB); % starting level spacing
            while Pdiff > tol && n < maxIterations  
                [~, noise_std]  = ber_apd_awgn(mpam, Tx, Fiber, this, Rx, sim);
                % noise_std assumes that levels and decision thresholds are
                % referred to the receiver

                mpam = mpam.optimize_level_spacing_gauss_approx(sim.BERtarget, Tx.Mod.rexdB, noise_std); % Optimized levels are with respect to transmitter
                % Levels optimized at the receiver
                
                % Required power at the APD input
                Pmean(n+1) = mean(mpam.a)/this.Geff;
                Tx.Ptx = Pmean(n+1);
                Pdiff = abs((Pmean(n+1) - Pmean(n))/Pmean(n));
                n = n+1; 
            end
            
            Pmean = Pmean(end);
            
            if n >= maxIterations
                warning('apd/optimize_PAM_levels: optimization did not converge')
            end
        end                
    end
           
    methods (Access=private)     
        function ber = calc_apd_ber(this, PtxdBm, Gapd, mpam, Tx, Fiber, Rx, sim)
            %% Iterate BER calculation: for given PtxdBm and Gapd calculates BER
            % This function is only called when M-PAM has equally spaced
            % levels
            Tx.Ptx = dBm2Watt(PtxdBm);

            % Set APD gain
            this.Gain = Gapd; % linear units

            ber = ber_apd_enumeration(mpam, Tx, Fiber, this, Rx, sim);
        end
    end
    
    methods
        %% Get and Set Methods
        function Fa = get.Fa(this) % excess noise factor
            %% Calculate Fa = Excess noise factor i.e., APD noise figure
            Fa = this.ka*this.Gain + (1 - this.ka)*(2 - 1/this.Gain); % Agrawal 4.4.18 (4th edition)
        end
               
        function GaindB = get.GaindB(this) 
            %% Calculate gain in dB
            GaindB = 10*log10(this.Gain);
        end
        
        function BW = get.BW(this)
            %% APD bandwidth
            BW = min(this.BW0, this.GainBW/this.Gain);
            % if this.GainBW/this.Gain < this.BW0, the APD is in the
            % avalanche buildup time limited regime
        end
        
        function Geff = get.Geff(this)
            %% Effective gain = Responsivity x Gain
            Geff = this.Gain*this.R;
        end
        
        function b = get.b(this) 
            %% Implicit relations of APD: beta = 1/(1 - ka)
            b = 1/(1-this.ka);
        end
        
        function a = get.a(this) 
            %% Implicit relations of APD: G = 1/(1 - ab)
            a =  1/this.b*(1-1/this.Gain);
        end
               
        %% Set methods
        function this = set.GaindB(this, GdB)
            %% Set gain in dB
            this.Gain = 10^(GdB/10); % set Gain, since GaindB is dependent
        end
        
        function self = setGain(self, Gain)
            %% Set APD Gain
            self.Gain = Gain;
        end
    end
end


%%%% Deprected methods %%%%
% Functions for modeling non-Gaussian noise statistics. These methods are
% no longer compatible with many functions they call.
%     %% Methods for calculating the accurate noise statistics (DEPRECTED)
%     methods
% %         function M = Ms(this, s)
% %             options = optimoptions('fsolve', 'Display', 'off');
% %             exitflag = 2;
% %             k = 1;
% %             while exitflag ~= 1 && k < 10
% %                 [M, ~, exitflag] = fsolve(@(M) M*(1 + this.a*(M-1))^(-this.b) - exp(s), 4*randn(1), options);
% %                 k = k + 1;
% %             end
% %             if exitflag ~= 1 
% %                 warning('Calculation of M did not converge');
% %             end  
% %         end       
% 
%         function M = Ms(this, s)
%             %% 
%             options = optimoptions('fsolve', 'Display', 'off');
%             for k = 1:length(s)
%                 [M(k), ~, exitflag] = fsolve(@(M) M*(1 + this.a*(M-1))^(-this.b) - exp(s(k)), 2*sign(s(k)), options);
%                 if exitflag ~= 1 
%                     warning('Calculation of M(s) did not converge');
%                 end
%             end  
%         end       
%         
%         %% Tail probabilities calculation
%         % Not working properly. fsolve doesn't work as well as fzero
%         function [px, shat] = tail_saddlepoint_approx(this, xthresh, P, fs, N0, tail) 
%             %% Output sgnal distribution including thermal noise using the saddlepoint approximation
%             % px = output sgnal pdf
%             % x = current at output
%             % lambda = rate of the Poisson process
%             % N0 = thermal noise psd
%             % fs = sampling frequency
%             %% NOT FINISHED !!!
%             options = optimoptions('fsolve', 'Display', 'off');
%             
%             % From implicit relations of APD: beta = 1/(1 - ka) and G = 1/(1 - ab)
%             a = this.a;
%             b = this.b;
% 
%             if P == 0
%                 px = 1;
%                 return;
%             end
%             
%             lambda = this.lambda(P, 1/fs);
% 
% %             if strcmp(tail, 'right')
% %                 sgn = 1;
% %             elseif strcmp(tail, 'left')
% %                 sgn = -1;
% %             else
% %                 error('invalid option');
% %             end
%            
%             dt = 1/fs; % pulse width
%             
%             % Calculate electron number variance corresponding to thermal noise
%             varTherm = (N0*fs/2)*(dt/this.q)^2;
%             
%             %
%             nthresh = xthresh*(dt/this.q);
% 
%             % dM(e^s)/ds
%             dMds = @(M) M*(1 + a*(M-1))/((1 + a*(M-1) - a*b*M)); 
% 
%             %% Calculate Saddlepoint
%             % Solve value of s 
%             % Can't use fzero because objective function might be complex
%             % Use real(s) instead of s to force real solution
%             % Solve value of M(e^s) at the saddlepoint 
%             % Can't use fzero because M(e^s) might be complex
%             [Mhat, ~, exitflag] = fsolve(@(M) lambda*(M*(1 + a*(M-1))/(1 + a*(M-1) - a*b*M))...
%                 + varTherm*(log(M) - b*log(1 + a*(M-1))) - nthresh - 1/(log(M) - b*log(1 + a*(M-1))), 1e-2, options);
% 
%             if exitflag ~= 1 % if didn't converge try again with opposite sgn
%                 [Mhat, ~, exitflag] = fsolve(@(M) lambda*(M*(1 + a*(M-1))/(1 + a*(M-1) - a*b*M))...
%                     + varTherm*(log(M) - b*log(1 + a*(M-1))) - nthresh - 1/(log(M) - b*log(1 + a*(M-1))), -1e-2, options);
%                 if exitflag ~= 1
%                     warning('Calculation of M at the saddlepoint did not converge');
%                 end
%             end
%             
%             % From M(e^s) at the saddlepoint calculate real part of s
%             shat = real(log(Mhat) - b*log(1 + a*(Mhat-1)));
%             
%             % Get new M(shat)
%             [Mhat, ~, exitflag] = fsolve(@(M) M*(1 + a*(M-1))^(-b) - exp(shat), real(Mhat), options);
%             if exitflag ~= 1 
%                 warning('Calculation of M at the saddlepoint did not converge');
%             end          
% 
%             
%             % First derivative of M(s) at the saddlepoint
%             dMhat = dMds(Mhat);
% 
%             % Second derivative of M(e^(s)) evaluated at the saddlepoint
%             d2Mhatds = d2Mds2(Mhat, dMhat);
% 
%             % 
%             Ksp = lambda*(Mhat - 1) + 1/2*varTherm*shat^2 - shat*nthresh - log(abs(shat));
%             d2Ksp = lambda*d2Mhatds + varTherm + 1/shat^2; % second derivative
% 
%             % Saddle point approximation
%             px = real(exp(Ksp)./sqrt(2*pi*d2Ksp));  
%                        
%             % Second derivative of M(e^s)
%             % M = M(e^s)
%             % dM = dM(e^s)/ds
%             function d2M = d2Mds2(M, dM)
%                 % M'' = (p'q - q'p)/q^2
% 
%                 p = M.*(1 + a*(M-1));
%                 pprime = dM.*(1 + a*(2*M-1));
%                 qq = (1 + a*(M-1) - a*b*M);
%                 qprime = a*dM*(1 - b);
% 
%                 % Second derivative of M(e^s)
%                 d2M = (pprime.*qq - qprime.*p)./qq.^2;
%             end
%         end
% 
%         function [px, x] = output_pdf_saddlepoint(this, P, fs, N0)
%             %% Output sgnal distribution including thermal noise using the saddlepoint approximation
%             % px = output sgnal pdf
%             % x = current at output
%             % lambda = rate of the Poisson process
%             % N0 = thermal noise psd
%             % fs = sampling frequency
%             options = optimoptions('fsolve', 'Display', 'off');
%             
%             % From implicit relations of APD: beta = 1/(1 - ka) and G = 1/(1 - ab)
%             a = this.a;
%             b = this.b;
% 
%             if P == 0
%                 px = 1;
%                 x = 0;
%                 return;
%             end
%             
%             % Rate of Poisson process
%             lambda = this.lambda(P, 1/fs);
% 
%             dt = 1/fs; % pulse width
%             
%             % Calculate electron number variance corresponding to thermal noise
%             varTherm = (N0*fs/2)*(dt/this.q)^2;
% 
%             % dM(e^s)/ds
%             dMds = @(M) M*(1 + a*(M-1))/((1 + a*(M-1) - a*b*M));
% 
%             % Uses gaussian appproximation to estimate number of points sufficient to obtain great part of the pdf
%             % this.Ptail is the probability of the tails not spanned by the
%             % calculated pdf
%             nmean = lambda*this.Gain;
%             nvar = this.Gain^2*this.Fa*lambda + varTherm;
%             npos = ceil(sqrt(nvar)*qfuncinv(this.Ptail));
%             N = max(0, nmean-npos):(nmean+npos);
% 
%             %% Calculate Saddlepoint
%             px = zeros(size(N));
%             for k = 1:length(N)
%                 n = N(k);
%                 
%                % Solve value of M(e^s) at the saddlepoint 
%                 % Can't use fzero because M(e^s) might be complex
%                 [Mhat, ~, exitflag] = fsolve(@(M) lambda*(M*(1 + a*(M-1))/(1 + a*(M-1) - a*b*M)) + varTherm*(log(M) - b*log(1 + a*(M-1))) - n, 1, options);
% 
%                 if exitflag ~= 1 % if didn't converge try again with opposite sgn
%                     [Mhat, ~, exitflag] = fsolve(@(M) lambda*(M*(1 + a*(M-1))/(1 + a*(M-1) - a*b*M)) + varTherm*(log(M) - b*log(1 + a*(M-1))) - n, -1, options);
%                     if exitflag ~= 1
%                         warning('Calculation of M at the saddlepoint did not converge');
%                         continue
%                     end
%                 end
%                
%                 % From M(e^s) at the saddlepoint calculate real part of s
%                 shat = real(log(Mhat) - b*log(1 + a*(Mhat-1)));
%                 
% %                 Mhat = this.Ms(shat);
% 
%                 % Get new M(shat)
%                 [Mhat, ~, exitflag] = fsolve(@(M) M*(1 + a*(M-1))^(-b) - exp(shat), real(Mhat), options);
%                 if exitflag ~= 1 
%                     warning('Calculation of M at the saddlepoint did not converge');
%                 end  
%                 
%                 % First derivative of M(s) at the saddlepoint
%                 dMhat = dMds(Mhat);
% 
%                 % Second derivative of M(e^(s)) evaluated at the saddlepoint
%                 d2Mhatds = d2Mds2(Mhat, dMhat);
% 
%                 % 
%                 Ksp = lambda*(Mhat - 1) + 1/2*varTherm*shat^2 - shat*n;
%                 d2Ksp = lambda*d2Mhatds + varTherm; % second derivative
% 
%                 % Saddle point approximation
%                 px(k) = real(exp(Ksp)./sqrt(2*pi*d2Ksp));  
%             end
% 
%             x = N*(this.q/dt);
%             px = px*(dt/this.q);
%             
%             intpx = trapz(x, px);
%             
%             if intpx < 1 - 2*this.Ptail
%                 warning(sprintf('output_pdf_saddlepoint: pdf accounts only for %g of probability, %g was expected', intpx, 1 - 2*this.Ptail));
%             end
%             
%             % renormalized pdf
%             px = px/intpx;
%                                     
%             % Second derivative of M(e^s)
%             % M = M(e^s)
%             % dM = dM(e^s)/ds
%             function d2M = d2Mds2(M, dM)
%                 % M'' = (p'q - q'p)/q^2
% 
%                 p = M.*(1 + a*(M-1));
%                 pprime = dM.*(1 + a*(2*M-1));
%                 qq = (1 + a*(M-1) - a*b*M);
%                 qprime = a*dM*(1 - b);
% 
%                 % Second derivative of M(e^s)
%                 d2M = (pprime.*qq - qprime.*p)./qq.^2;
%             end            
%         end
%         
%         % Calculate noise pdf for a sgnal levels Plevels with duration dt.
%         % The Gaussian approximation is compared with the distribution 
%         % calculated using the saddlepoint approximation
%         function lpdf = levels_pdf(this, Plevels, fs)           
%             % Struct of levels pdf: 
%             lpdf = struct('p', 0, 'p_gauss', 0, 'I', 0,... % p(I) = true pdf, p_gauss(I) = pdf under Gaussian approximation
%                 'mean', 0, 'mean_gauss', 0,... % true mean and mean under Gaussian approximation
%                 'var', 0, 'var_gauss', num2cell(zeros(size(Plevels)))); % true variance and variance under Gaussian approximation
%             
%             for k = 1:length(Plevels)                  
%                 
%                 if Plevels(k) == 0
%                     continue
%                 end
%                 
%                 [lpdf(k).p, lpdf(k).I] = this.output_pdf_saddlepoint(Plevels(k), fs, 0);
%                                                
%                 lpdf(k).mean = trapz(lpdf(k).I, lpdf(k).I.*lpdf(k).p);
%                 lpdf(k).mean_gauss = Plevels(k)*this.Gain;
%                 lpdf(k).var = trapz(lpdf(k).I, lpdf(k).I.^2.*lpdf(k).p) - lpdf(k).mean.^2;
%                 lpdf(k).var_gauss = this.varShot(Plevels(k), fs/2);
%                 
%                 lpdf(k).p_gauss = pdf('normal', lpdf(k).I, lpdf(k).mean_gauss, sqrt(lpdf(k).var_gauss));
%             end
%         end
%         
%         %% Output sgnal pmf (without thermal noise)
%         % pn = probability of observing n electrons at the output
%         % lambda = rate of the Poisson process
%         function pn = output_pmf(this, lambda)          
%             pn = exp(-lambda); % n = 0;
%             
%             psum = pn;
%             k = 0; % current iteration   
%             r = 0;
%             mr1 = this.calc_mr(r+1);
%             while psum  < this.cdf_accuracy && k < this.Niterations              
%                 pn(k+2) = lambda/(k+1)*sum((r+1).*mr1.*fliplr(pn(r+1)));
%                 
%                 psum = psum + pn(k+2);
%                 
%                 k = k + 1;
%                 
%                 r = [r k];
%                 
%                 mr1 = [mr1 calc_mr(k+1)];
%             end
%             
%             if k >= this.Niterations
%                 warning(sprintf('output_pmf(this, lambda): max number of iterations exceeded. pmf accounts only for %f of probability', psum));
%             end
%             
%             % Calculate mr given in C. Helstrom "Computattion of Output
%             % Electron Distributions in Avalanche Photodiodes"
%             function mr = calc_mr(r)
% 
%                 P2 = (r-1)*log(this.a);
%                 P3 = (r*(this.b-1)+1)*log(1 - this.a);
%                 P4 = sum(log(1:r)); % P4 = log(factorial(r));
% 
%                 if r < 100 % calculate exactly
%                     P1 = log(gamma(this.b*r + 1));
%                     P5 = log(gamma(r*(this.b-1) + 2));
%                 else % use Stirling's formula for the factorial and Gamma functions
%                     stirling_approx = @(z) log(sqrt(2*pi/z)) + z*log(z/exp(1));
%                     P1 = stirling_approx(this.b*r + 1);
%                     P5 = stirling_approx(r*(this.b-1) + 2);                
%                 end      
% 
%                 P = P1 + P2 + P3 - P4 - P5;
% 
%                 mr = exp(P);
%             end
%         end  
%     end