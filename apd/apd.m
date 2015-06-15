%% Avalanche photodiode
classdef apd < handle
    properties
        Gain % Gain
        ka   % impact ionization factor
        GainBW % Gain x Bandwidth product
        R    % responsivity
        Id   % dark current
    end
    properties (Dependent)
        Fa % excess noise factor
        GaindB
    end
    
    properties (Constant)
        h = 6.62606957e-34; % Planck
        q = 1.60217657e-19; % electron charge
        c = 299792458;      % speed of light
    end
    
    properties (Constant, GetAccess=private)
        cdf_accuracy = 1-1e-4; % Required accuracy of the cdf (i.e., pmf will have the minimum number of points that satistifes that sum pmf > cdf_accuracy.)
        Niterations = 1e6; % maximum number of iterations in a while loop
        Ptail = 1e-6; % probability of clipped tail
    end
    
    properties (Dependent, GetAccess=private)
        a % auxiliary variable G = 1/(1-ab)
        b % auxiliary variable b = 1/(1-ka)
    end
    
    methods
        %% constructor
        function this = apd(GaindB, ka, GainBW, R, Id) 
            this.GainBW = GainBW;
            this.Gain = 10^(GaindB/10);
            this.ka = ka;
                        
            if nargin == 4
                this.R = R;
            else 
                this.R = 1;
            end
            
            if nargin == 5
                this.Id = Id;
            else 
                this.Id = 0;
            end
           
        end
                           
        %% Excess noise factor
        function Fa = get.Fa(this) % Excess noise
            Fa = this.ka*this.Gain + (1 - this.ka)*(2 - 1/this.Gain); % Agrawal 4.4.18 (4th edition)
        end
               
        function GaindB = get.GaindB(this) % Excess noise
            GaindB = 10*log10(this.Gain);
        end
        
        function b = get.b(this)
            b = 1/(1-this.ka);
        end
        
        function a = get.a(this)
            a =  1/this.b*(1-1/this.Gain);
        end
                     
        %% Detection
        % Pin = received power; fs = sampling frequency
        % noise_stats = 'gaussian' or 'doubly-stochastic'
        % N0 = thermal noise psd
        function output = detect(this, Pin, fs, noise_stats, N0)
            if nargin < 5
                N0 = 0;
            end                
            
            if strcmp(noise_stats, 'gaussian')
                % Assuming Gaussian statistics
                var_shot = 2*this.q*this.Gain^2*this.Fa*(this.R*Pin + this.Id)*fs/2; % Agrawal 4.4.17 (4th edition)
                
                var_therm = N0*fs/2; % thermal noise
                
                output = this.R*this.Gain*Pin + sqrt(var_shot + var_therm).*randn(size(Pin));
              
            % uses saddlepoint approximation to obtain pmf in order to generate 
            % output distributed according to that pmf
            elseif  strcmp(noise_stats, 'doubly-stochastic')  
                Plevels = unique(Pin);
                
                % pulse width
                dt = 1/fs; 
               
                output = zeros(size(Pin));
                for k = 1:length(Plevels)
                    lambda = (this.R*Plevels(k) + this.Id)*dt/this.q; % Personick, "Statistics of a General Class of Avalanche Detectors With Applications to Optical Communication"
                    [px, x] = this.output_pdf_saddlepoint(lambda, fs, 0); % doesn't include thermal noise here
                    
                    cdf = cumtrapz(x, px);
                                                   
                    pos = (Pin == Plevels(k));
                    
                    % Sample according to pmf pk_P
                    u = rand(sum(pos), 1); % uniformly-distributed
                    dist = abs(bsxfun(@minus, u, cdf));
                    [~, ix] = min(dist, [], 2);
                    output(pos) = x(ix); % distributed accordingly to px

                end
                
                output = output + sqrt(N0*fs/2).*randn(size(Pin)); % includes thermal noise
            else
                error('Invalid Option!')
            end
        end
        
        % Optimize apd gain for a given system
        function optimize_gain(this, mpam, tx, rx, sim)

            if strcmp(mpam.level_spacing, 'uniform')
                % Optmize gain for uniform spacing: find Gapd that minimizes the required
                % average power (Prec) to achieve a certain target SER.
                [Gapd_opt, ~, exitflag] = fminbnd(@(Gapd) fzero(@(PtxdBm) calc_apd_ber(PtxdBm, Gapd, mpam, tx, this, rx, sim) - sim.BERtarget, -20), 1, min(this.GainBW/mpam.Rs, 100));    

            elseif strcmp(mpam.level_spacing, 'nonuniform')
                % Optimal level spacing
                [Gapd_opt, ~, exitflag] = fminbnd(@(Gapd) mean(calc_level_spacing(Gapd, mpam, tx, this, rx, sim))/Gapd, 1, min(this.GainBW/mpam.Rs, 100));
                1;
            else
                error('Invalid Option')
            end

            if exitflag ~= 1
                warning(sprintf('APD gain optimization did not converge (exitflag = %d)\n', exitflag))
            end 

            if ~isnan(Gapd_opt) && ~isinf(Gapd_opt)
                this.Gain = Gapd_opt;
            else
                warning('apd_gain_optmization: APD gain was not changed')
            end

            function ber = calc_apd_ber(PtxdBm, Gapd, mpam, tx, apd, rx, sim)
                % Set power level
                tx.PtxdBm = PtxdBm;

                % Set APD gain
                apd.Gain = Gapd; % linear units

                ber_struct = apd_ber(mpam, tx, apd, rx, sim);

                ber = ber_struct.gauss;
            end

            function a = calc_level_spacing(Gapd, mpam, tx, apd, rx, sim)
                apd.Gain = Gapd; % linear units

                a = level_spacing_optm_gauss_approx(mpam, tx, apd, rx, sim);
            end

        end

        %% Output sgnal pmf (without thermal noise)
        % pn = probability of observing n electrons at the output
        % lambda = rate of the Poisson process
        function pn = output_pmf(this, lambda)          
            pn = exp(-lambda); % n = 0;
            
            psum = pn;
            k = 0; % current iteration   
            r = 0;
            mr1 = this.calc_mr(r+1);
            while psum  < this.cdf_accuracy && k < this.Niterations              
                pn(k+2) = lambda/(k+1)*sum((r+1).*mr1.*fliplr(pn(r+1)));
                
                psum = psum + pn(k+2);
                
                k = k + 1;
                
                r = [r k];
                
                mr1 = [mr1 calc_mr(k+1)];
            end
            
            if k >= this.Niterations
                warning(sprintf('output_pmf(this, lambda): max number of iterations exceeded. pmf accounts only for %f of probability', psum));
            end
            
            % Calculate mr given in C. Helstrom "Computattion of Output
            % Electron Distributions in Avalanche Photodiodes"
            function mr = calc_mr(r)

                P2 = (r-1)*log(this.a);
                P3 = (r*(this.b-1)+1)*log(1 - this.a);
                P4 = sum(log(1:r)); % P4 = log(factorial(r));

                if r < 100 % calculate exactly
                    P1 = log(gamma(this.b*r + 1));
                    P5 = log(gamma(r*(this.b-1) + 2));
                else % use Stirling's formula for the factorial and Gamma functions
                    stirling_approx = @(z) log(sqrt(2*pi/z)) + z*log(z/exp(1));
                    P1 = stirling_approx(this.b*r + 1);
                    P5 = stirling_approx(r*(this.b-1) + 2);                
                end      

                P = P1 + P2 + P3 - P4 - P5;

                mr = exp(P);
            end
        end
                
        %% Output sgnal distribution including thermal noise using the saddlepoint approximation
        % px = output sgnal pdf
        % x = current at output
        % lambda = rate of the Poisson process
        % N0 = thermal noise psd
        % fs = sampling frequency
        function [px, shat] = tail_saddlepoint_approx(this, xthresh, lambda, fs, N0, tail)
            options = optimoptions('fsolve', 'Display', 'off');
            
            % From implicit relations of APD: beta = 1/(1 - ka) and G = 1/(1 - ab)
            a = this.a;
            b = this.b;

            if lambda == 0
                px = 1;
                return;
            end

            if strcmp(tail, 'right')
                sgn = 1;
            elseif strcmp(tail, 'left')
                sgn = -1;
            else
                error('invalid option');
            end
           
            dt = 1/fs; % pulse width
            
            % Calculate electron number variance corresponding to thermal noise
            varTherm = (N0*fs/2)*(dt/this.q)^2;
            
            %
            nthresh = xthresh*(dt/this.q);

            % dM(e^s)/ds
            dMds = @(M) M*(1 + a*(M-1))/((1 + a*(M-1) - a*b*M)); 

            %% Calculate Saddlepoint
            shat = -sgn;
            exitflag = 2;
            while sign(shat) ~= sgn && exitflag ~= 1
                % Solve value of M(e^s) at the saddlepoint 
                % Can't use fzero because M(e^s) might be complex
                [Mhat, 	~, exitflag] = fsolve(@(M) lambda*(M*(1 + a*(M-1))/(1 + a*(M-1) - a*b*M))...
                    + varTherm*(log(M) - b*log(1 + a*(M-1))) - nthresh - 1/(log(M) - b*log(1 + a*(M-1))), 2*randn(1), options);

                if exitflag ~= 1
                    warning('tail_saddlepoint_approx: Did not converge to a solution of M (1)');
                end

                % From M(e^s) at the saddlepoint calculate real part of s
                shat = real(log(Mhat) - b*log(1 + a*(Mhat-1)));
            end
            
%             if sign(shat) ~= sgn
%                 % Solve value of M(e^s) at the saddlepoint 
%                 % Can't use fzero because M(e^s) might be complex
%                 [Mhat, fval, exitflag] = fsolve(@(M) lambda*(M*(1 + a*(M-1))/(1 + a*(M-1) - a*b*M))...
%                     + varTherm*(log(M) - b*log(1 + a*(M-1))) - nthresh - 1/(log(M) - b*log(1 + a*(M-1))), sgn*2, options);
% 
%                 if exitflag ~= 1 % if didn't converge try again with opposite sgn
%                     [Mhat, fval, exitflag] = fsolve(@(M) lambda*(M*(1 + a*(M-1))/(1 + a*(M-1) - a*b*M))...
%                         + varTherm*(log(M) - b*log(1 + a*(M-1))) - nthresh - 1/(log(M) - b*log(1 + a*(M-1))), -sgn*2, options);
% 
%                     if exitflag ~= 1
%                         warning('tail_saddlepoint_approx: Did not converge to a solution of M (2)');
%                     end
%                 end
%             end
% 
%             % From M(e^s) at the saddlepoint calculate real part of s
%             shat = real(log(Mhat) - b*log(1 + a*(Mhat-1)));
        
            % Get new M(shat)
            [Mhat, ~, exitflag] = fsolve(@(M) M*(1 + a*(M-1))^(-b) - exp(shat), real(Mhat), options);
            if exitflag ~= 1 
                warning('Calculation of M at the saddlepoint did not converge');
            end            

            % First derivative of M(s) at the saddlepoint
            dMhat = dMds(Mhat);

            % Second derivative of M(e^(s)) evaluated at the saddlepoint
            d2Mhatds = d2Mds2(Mhat, dMhat);

            % 
            Ksp = lambda*(Mhat - 1) + 1/2*varTherm*shat^2 - shat*nthresh - log(abs(shat));
            d2Ksp = lambda*d2Mhatds + varTherm + 1/shat^2; % second derivative

            % Saddle point approximation
            px = real(exp(Ksp)./sqrt(2*pi*d2Ksp));  
                       
            % Second derivative of M(e^s)
            % M = M(e^s)
            % dM = dM(e^s)/ds
            function d2M = d2Mds2(M, dM)
                % M'' = (p'q - q'p)/q^2

                p = M.*(1 + a*(M-1));
                pprime = dM.*(1 + a*(2*M-1));
                qq = (1 + a*(M-1) - a*b*M);
                qprime = a*dM*(1 - b);

                % Second derivative of M(e^s)
                d2M = (pprime.*qq - qprime.*p)./qq.^2;
            end
        end

        %% Output sgnal distribution including thermal noise using the saddlepoint approximation
        % px = output sgnal pdf
        % x = current at output
        % lambda = rate of the Poisson process
        % N0 = thermal noise psd
        % fs = sampling frequency
        function [px, x] = output_pdf_saddlepoint(this, lambda, fs, N0)
            options = optimoptions('fsolve', 'Display', 'off');
            
            % From implicit relations of APD: beta = 1/(1 - ka) and G = 1/(1 - ab)
            a = this.a;
            b = this.b;

            if lambda == 0
                px = 1;
                x = 0;
                return;
            end

            dt = 1/fs; % pulse width
            
            % Calculate electron number variance corresponding to thermal noise
            varTherm = (N0*fs/2)*(dt/this.q)^2;

            % 
            dMds = @(M) M*(1 + a*(M-1))/((1 + a*(M-1) - a*b*M)); %  dM(e^s)/ds

            % Uses gaussian appproximation to estimate number of points sufficient to obtain great part of the pdf
            nmean = lambda*this.Gain;
            nvar = this.Gain^2*this.Fa*lambda + varTherm;
            npos = ceil(sqrt(nvar)*qfuncinv(this.Ptail));
            N = max(0, nmean-npos):(nmean+npos);

            %% Calculate Saddlepoint
            px = zeros(size(N));
            for k = 1:length(N)
                n = N(k);

                % Solve value of M(e^s) at the saddlepoint 
                % Can't use fzero because M(e^s) might be complex
                [Mhat, fval, exitflag] = fsolve(@(M) lambda*(M*(1 + a*(M-1))/(1 + a*(M-1) - a*b*M)) + varTherm*(log(M) - b*log(1 + a*(M-1))) - n, 1, options);

                if exitflag ~= 1 % if didn't converge try again with opposite sgn
                    [Mhat, fval, exitflag] = fsolve(@(M) lambda*(M*(1 + a*(M-1))/(1 + a*(M-1) - a*b*M)) + varTherm*(log(M) - b*log(1 + a*(M-1))) - n, -1, options);
                    if exitflag ~= 1
                        warning('Calculation of M at the saddlepoint did not converge');
                        continue
                    end
                end
               
                % From M(e^s) at the saddlepoint calculate real part of s
                shat = real(log(Mhat) - b*log(1 + a*(Mhat-1)));
                               
                % Get new M(shat)
                [Mhat, ~, exitflag] = fsolve(@(M) M*(1 + a*(M-1))^(-b) - exp(shat), real(Mhat), options);
                if exitflag ~= 1 
                    warning('Calculation of M at the saddlepoint did not converge');
                end            

                % First derivative of M(s) at the saddlepoint
                dMhat = dMds(Mhat);

                % Second derivative of M(e^(s)) evaluated at the saddlepoint
                d2Mhatds = d2Mds2(Mhat, dMhat);

                % 
                Ksp = lambda*(Mhat - 1) + 1/2*varTherm*shat^2 - shat*n;
                d2Ksp = lambda*d2Mhatds + varTherm; % second derivative

                % Saddle point approximation
                px(k) = real(exp(Ksp)./sqrt(2*pi*d2Ksp));  
            end

            x = N*(this.q/dt);
            px = px*(dt/this.q);
            
            intpx = trapz(x, px);
            
            if intpx < 1 - 2*this.Ptail
                warning(sprintf('output_pdf_saddlepoint: pdf accounts only for %g of probability, %g was expected', intpx, 1 - 2*this.Ptail));
            end
            
            % renormalized pdf
            px = px/intpx;
                                    
            % Second derivative of M(e^s)
            % M = M(e^s)
            % dM = dM(e^s)/ds
            function d2M = d2Mds2(M, dM)
                % M'' = (p'q - q'p)/q^2

                p = M.*(1 + a*(M-1));
                pprime = dM.*(1 + a*(2*M-1));
                qq = (1 + a*(M-1) - a*b*M);
                qprime = a*dM*(1 - b);

                % Second derivative of M(e^s)
                d2M = (pprime.*qq - qprime.*p)./qq.^2;
            end
        end
        
        % Calculate noise pdf for a sgnal levels Plevels with duration dt.
        % The Gaussian approximation is compared with the distribution 
        % calculated using the saddlepoint approximation
        function lpdf = levels_pdf(this, Plevels, fs)
            dt = 1/fs; % pulse width
            
            % Struct of levels pdf: 
            lpdf = struct('p', 0, 'p_gauss', 0, 'I', 0,... % p(I) = true pdf, p_gauss(I) = pdf under Gaussian approximation
                'mean', 0, 'mean_gauss', 0,... % true mean and mean under Gaussian approximation
                'var', 0, 'var_gauss', num2cell(zeros(size(Plevels)))); % true variance and variance under Gaussian approximation
            
            for k = 1:length(Plevels)                  
                
                if Plevels(k) == 0
                    continue
                end
                
                lambda = (this.R*Plevels(k) + this.Id)*dt/this.q; % Personick, "Statistics of a General Class of Avalanche Detectors With Applications to Optical Communication"
                
                [lpdf(k).p, lpdf(k).I] = this.output_pdf_saddlepoint(lambda, fs, 0);
                                               
                lpdf(k).mean = trapz(lpdf(k).I, lpdf(k).I.*lpdf(k).p);
                lpdf(k).mean_gauss = Plevels(k)*this.Gain;
                lpdf(k).var = trapz(lpdf(k).I, lpdf(k).I.^2.*lpdf(k).p) - lpdf(k).mean.^2;
                lpdf(k).var_gauss = 2*this.q*this.Gain^2*this.Fa*(this.R*Plevels(k) + this.Id)*fs/2; 
                
                lpdf(k).p_gauss = pdf('normal', lpdf(k).I, lpdf(k).mean_gauss, sqrt(lpdf(k).var_gauss));
            end
        end
        
    end
end
    


% 
%   %% Random Gain distribution using the saddle point approximation
%         %% Does not include thermal noise
%         % pn = probability of observing n electrons at the output
%         % n = number of electrons at the output
%         % lambda = rate of the Poisson process
%         % Does not include thermal noise
%         function [pn, n] = output_pmf_saddlepoint(this, lambda)
%             a = this.a;
%             b = this.b;
%             
%             if lambda == 0
%                 pn = 1;
%                 n = 0;
%                 return;
%             end
%             
%             % Estimate number of points sufficient to obtain great part of the pmf
%             nmean = lambda*this.Gain;
%             nvar = this.Gain^2*this.Fa*lambda;
%             npos = fzero(@(x) qfunc(x/sqrt(nvar)) - this.Ptail, 5*sqrt(nvar));
%             n = max(0, nmean-npos):(nmean+npos);
% %             n = 0:(lambda*this.Gain*5);
%             
%             %% Calculate Saddle-point
%             % based on the derivation by C. W. Helstrom in "Computation of
%             % Output Electron Distributions in Avalanche Photodiodes"
%             tau = n + 1;
%             gamma = lambda*(1-a) + a*(b - 1)*tau;
%             
%             % probability generating function of p(n|1) evaluated at the
%             % saddle point
%             Msp = (sqrt(gamma.^2 + 4*a*(1-a)*lambda*tau) - gamma)/(2*a*lambda);
%             
%             % z at the saddle point (should be real)
%             zsp = Msp.*(1 + a*(Msp-1)).^(-b);
%             
%             % First and second derivative of M(z) evaluated at the saddle point
%             % M' = p/q
%             % M'' = (p'q - q'p)/q^2
%             
%             M = Msp; 
%             dM = M.*(1 + a*(M-1))./(zsp.*(1 + a*(M - 1) - a*b*M));
%             
%             p = M.*(1 + a*(M-1));
%             pprime = dM.*(1 + a*(2*M-1));
%             qq = zsp.*(1 + a*(M-1) - a*b*M);
%             qprime = 1 + a*(M + zsp.*dM - 1) - a*b*(M + zsp.*dM);
%             
%             % Second derivative of M(z) evaluated at the saddle point
%             ddM = (pprime.*qq - qprime.*p)./qq.^2;
%             
%             % Phase function (i.e., Ksp = lambda(M(z)-1) - nln(z)
%             Ksp = lambda*(M -1) - n.*log(zsp);
%             ddKsp = lambda*ddM + n./zsp.^2; % second derivative
%             
%             % Saddle point approximation
%             pn = exp(Ksp)./sqrt(2*pi*ddKsp);  
%         end




            