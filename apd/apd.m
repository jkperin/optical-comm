%% Avalanche photodiode
classdef apd
    properties
        GainBW % Gain x Bandwidth product
        Gain % Gain
        ka   % impact ionization factor
        R    % responsivity
        Id   % dark current
        noise_stats % shot noise statistics: {'off', 'Gaussian', 'DoublyStoch'}, default = Gaussian
    end
    properties (Dependent)
        Fa % excess noise factor
    end
    
    properties (Constant)
        h = 6.62606957e-34; % Planck
        q = 1.60217657e-19; % electron charge
        c = 299792458;      % speed of light
    end
    
    properties (Constant, GetAccess=private)
        cdf_accuracy = 1-1e-5; % Required accuracy of the cdf (i.e., pmf will have the minimum number of points that satistifes that sum pmf > cdf_accuracy.)
        Niterations = 1e6; % maximum number of iterations in a while loop
        Ptail = 1e-12; % probability of clipped tail
    end
    
    properties (Dependent, GetAccess=private)
        a % auxiliary variable G = 1/(1-ab)
        b % auxiliary variable b = 1/(1-ka)
    end
    
    methods
        %% constructor
        function obj = apd(GainBW, ka, Gain, noise_stats, R, Id) 
            obj.GainBW = GainBW;
            obj.Gain = Gain;
            obj.ka = ka;
            
            if nargin == 4
                obj.noise_stats = noise_stats;
            else
                obj.noise_stats = 'Gaussian';
            end
            
            if nargin == 5
                obj.R = R;
            else 
                obj.R = 1;
            end
            
            if nargin == 6
                obj.Id = Id;
            else 
                obj.Id = 0;
            end
           
        end
        
        %% Set noise statistics
        function obj = set.noise_stats(obj, str)
            if strcmp(str, 'Gaussian') || strcmp(str, 'DoublyStoch') || strcmp(str, 'off')
                obj.noise_stats = str;
            else
                error('Invalid option: noise_stats should be either off, Gaussian or DoublyStoch');
            end
        end
                    
        %% Excess noise factor
        function Fa = get.Fa(obj) % Excess noise
            Fa = obj.ka*obj.Gain + (1 - obj.ka)*(2 - 1/obj.Gain); % Agrawal 4.4.18 (4th edition)
        end
        
        function b = get.b(obj)
            b = 1/(1-obj.ka);
        end
        
        function a = get.a(obj)
            a =  1/obj.b*(1-1/obj.Gain);
        end
               
        %% Detection
        % Pin = received power; dt = sampling period: 1/dt = noise
        % bandwidth
        function output = detect(obj, Pin, dt)
            if strcmp(obj.noise_stats, 'Gaussian')
                % Assuming Gaussian statistics
                var_shot = 2*obj.q*obj.Gain^2*obj.Fa*(obj.R*Pin + obj.Id)*1/dt; % Agrawal 4.4.17 (4th edition)
                
                output = obj.R*obj.Gain*Pin + sqrt(var_shot).*randn(size(Pin));
            else 
                Plevels = unique(Pin);
                
                %% !!!
                dt = dt/2; % required to match gaussian approximation
               
                nelec = zeros(size(Pin));
                for k = 1:length(Plevels)
                    lambda = (obj.R*Plevels(k) + obj.Id)*dt/obj.q; % Personick, "Statistics of a General Class of Avalanche Detectors With Applications to Optical Communication"
                    [pk_P, n] = obj.output_pmf_saddlepoint(lambda);
                    
                    cdf = cumsum(pk_P);
                    
                    if cdf(end) <= obj.cdf_accuracy
                        warning(sprintf('output_pmf_saddlepoint(obj, lambda): pmf accounts only for %f of probability', cdf(end)));
                    end
                                
                    pos = (Pin == Plevels(k));
                    
                    % Sample according to pmf pk_P
                    u = rand(sum(pos), 1); % uniformly-distributed
                    dist = abs(bsxfun(@minus, u, cdf));
                    [~, ix] = min(dist, [], 2);
                    nelec(pos) = n(ix); % n distributed accordingly to p

                end
                
                output = nelec*obj.q/dt; 
            end
        end
       
        %% Random Gain distribution
        % pn = probability of observing n electrons at the output
        % lambda = rate of the Poisson process
        % lambda = rate of the Poisson process
        function pn = output_pmf(obj, lambda)          
            pn = exp(-lambda); % n = 0;
            
            psum = pn;
            k = 0; % current iteration   
            r = 0;
            mr1 = obj.calc_mr(r+1);
            while psum  < obj.cdf_accuracy && k < obj.Niterations              
                pn(k+2) = lambda/(k+1)*sum((r+1).*mr1.*fliplr(pn(r+1)));
                
                psum = psum + pn(k+2);
                
                k = k + 1;
                
                r = [r k];
                
                mr1 = [mr1 obj.calc_mr(k+1)];
            end
            
            if k >= obj.Niterations
                warning(sprintf('output_pmf(obj, lambda): max number of iterations exceeded. pmf accounts only for %f of probability', psum));
            end
            
        end
        
        %% Random Gain distribution using the saddle point approximation
        % pn = probability of observing n electrons at the output
        % n = number of electrons at the output
        % lambda = rate of the Poisson process
        function [pn, n] = output_pmf_saddlepoint(obj, lambda)
            a = obj.a;
            b = obj.b;
            
            if lambda == 0
                pn = 1;
                n = 0;
                return;
            end
            
            % Estimate number of points sufficient to obtain great part of the pmf
            nmean = lambda*obj.Gain;
            nvar = obj.Gain^2*obj.Fa*lambda;
            npos = fzero(@(x) qfunc(x/sqrt(nvar)) - obj.Ptail, 5*sqrt(nvar));
            n = max(0, nmean-npos):(nmean+npos);
%             n = 0:(lambda*obj.Gain*5);
            
            %% Calculate Saddle-point
            % based on the derivation by C. W. Helstrom in "Computation of
            % Output Electron Distributions in Avalanche Photodiodes"
            tau = n + 1;
            gamma = lambda*(1-a) + a*(b - 1)*tau;
            
            % probability generating function of p(n|1) evaluated at the
            % saddle point
            Msp = (sqrt(gamma.^2 + 4*a*(1-a)*lambda*tau) - gamma)/(2*a*lambda);
            
            % z at the saddle point (should be real)
            zsp = Msp.*(1 + a*(Msp-1)).^(-b);
            
            % First and second derivative of M(z) evaluated at the saddle point
            % M' = p/q
            % M'' = (p'q - q'p)/q^2
            
            M = Msp; 
            dM = M.*(1 + a*(M-1))./(zsp.*(1 + a*(M - 1) - a*b*M));
            
            p = M.*(1 + a*(M-1));
            pprime = dM.*(1 + a*(2*M-1));
            qq = zsp.*(1 + a*(M-1) - a*b*M);
            qprime = 1 + a*(M + zsp.*dM - 1) - a*b*(M + zsp.*dM);
            
            % Second derivative of M(z) evaluated at the saddle point
            ddM = (pprime.*qq - qprime.*p)./qq.^2;
            
            % Phase function (i.e., Ksp = lambda(M(z)-1) - nln(z)
            Ksp = lambda*(M -1) - n.*log(zsp);
            ddKsp = lambda*ddM + n./zsp.^2; % second derivative
            
            % Saddle point approximation
            pn = exp(Ksp)./sqrt(2*pi*ddKsp);  
        end
        
        
        % Calculate mr given in C. Helstrom "Computattion of Output
        % Electron Distributions in Avalanche Photodiodes"
        function mr = calc_mr(obj, r)
            
            P2 = (r-1)*log(obj.a);
            P3 = (r*(obj.b-1)+1)*log(1 - obj.a);
            P4 = sum(log(1:r)); % P4 = log(factorial(r));
            
            if r < 100 % calculate exactly
                P1 = log(gamma(obj.b*r + 1));
                P5 = log(gamma(r*(obj.b-1) + 2));
            else % use Stirling's formula for the factorial and Gamma functions
                stirling_approx = @(z) log(sqrt(2*pi/z)) + z*log(z/exp(1));
                P1 = stirling_approx(obj.b*r + 1);
                P5 = stirling_approx(r*(obj.b-1) + 2);                
            end      
            
            P = P1 + P2 + P3 - P4 - P5;
            
            mr = exp(P);
        end
        
        function lpdf = levels_pdf(obj, Plevels, dt)
            dt = dt/2; % required to match gaussian approximation
            
            for k = 1:length(Plevels)                  
                
                if Plevels == 0
                    lpdf(k).p = 0;
                    lpdf(k).I = 0;
                    lpdf(k).mean = 0;
                    lpdf(k).mean_gauss = 0;
                    lpdf(k).p_gauss = 0;
                    continue
                end
                
                lambda = (obj.R*Plevels(k) + obj.Id)*dt/obj.q; % Personick, "Statistics of a General Class of Avalanche Detectors With Applications to Optical Communication"
                
                [lpdf(k).p, lpdf(k).I] = obj.output_pmf_saddlepoint(lambda);
                
                lpdf(k).p = lpdf(k).p*(dt/obj.q);
                lpdf(k).I = lpdf(k).I*(obj.q/dt);
                               
                lpdf(k).mean = trapz(lpdf(k).I, lpdf(k).I.*lpdf(k).p);
                lpdf(k).mean_gauss = Plevels(k)*obj.Gain;
                lpdf(k).var = trapz(lpdf(k).I, lpdf(k).I.^2.*lpdf(k).p) - lpdf(k).mean.^2;
                lpdf(k).var_gauss = obj.q*obj.Gain^2*obj.Fa*(obj.R*Plevels(k) + obj.Id)*1/dt; % dt = dt/2 so factor of 2 is not necessary for one-sided psd
                
                lpdf(k).p_gauss = pdf('normal', lpdf(k).I, lpdf(k).mean_gauss, sqrt(lpdf(k).var_gauss));
            end
        end
        
    end
end
    

            