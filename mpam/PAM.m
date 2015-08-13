classdef PAM < handle
    properties
        M % constellation size
        Rb % bit rate
        level_spacing % 'equally-spaced' or 'optimized'
        pshape % pulse shape function
        a % levels
        b % decision threshold
    end
    
    properties (Dependent)
        Rs % symbol rate
    end
    
    properties (GetAccess=private)
        % Used in level spacing optimization
        maxtol = 1e-6; % maximum tolerance for convergence
        maxit = 20; % maximum number of iteratios
    end
       
    methods
        function obj = PAM(M, Rb, level_spacing, pshape)
            obj.M = M;
            obj.Rb = Rb;
            obj.level_spacing = level_spacing;
            obj.pshape = pshape;
            
            switch level_spacing
                case 'equally-spaced'
                    obj.a = ((0:2:2*(M-1))/(2*(M-1))).';
                    obj.b = ((1:2:(2*(M-1)-1))/(2*(M-1))).';
                case 'optimized'
                    % Optimize level spacing function must be called
                    obj.a = []; 
                    obj.b = [];
                otherwise
                    error('pam class: Invalid level spacing option')
            end
        end
        
        %% Get methods
        function Rs = get.Rs(this)
            Rs = this.Rb/log2(this.M);
        end
        
        % Set levels to desired values
        function set_levels(this, levels, thresholds)
            this.a = levels/levels(end);
            this.b = thresholds/levels(end);
        end
        
        % Normalize levels so that last level is unit
        function norm_levels(this)
            this.b = this.b/this.a(end);
            this.a = this.a/this.a(end);
        end
                    
        % Adjust levels to desired transmitted power and extinction ratio
        function [Plevels, Pthresh] = adjust_levels(this, Ptx, rexdB)
            rex = 10^(-abs(rexdB)/10); % extinction ratio. Defined as Pmin/Pmax
            switch this.level_spacing
                case 'equally-spaced'
                    % Restart levels
                    this.a = ((0:2:2*(this.M-1))/(2*(this.M-1))).';
                    this.b = ((1:2:(2*(this.M-1)-1))/(2*(this.M-1))).';
                    
                    amean = mean(this.a); % mean value
                    
                    Pmin = 2*Ptx*rex/(1 + rex); % power of the lowest level 
                    Plevels = this.a*(Ptx/amean)*((1-rex)/(1+rex)) + Pmin; % levels at the transmitter
                    Pthresh = this.b*(Ptx/amean)*((1-rex)/(1+rex)) + Pmin; % decision thresholds at the transmitter
                case 'optimized'
                    amean = mean(this.a); % mean value
                    
                    % Extinction ratio was already enforced in the
                    % optimization process, so just scale to desired power
                    Plevels = this.a*Ptx/amean; % levels at the transmitter
                    Pthresh = this.b*Ptx/amean; % decision thresholds at the transmitter
                otherwise
                    error('pam class: Invalid level spacing option')
            end
            
            this.a = Plevels;
            this.b = Pthresh;            
        end
        
        % Generate PAM signal
        function [xt, xd] = mod(this, dataTX, Mct)
            xd = this.a(gray2bin(dataTX, 'pam', this.M) + 1);
            xt = reshape(kron(xd, this.pshape(0:Mct-1)).', length(dataTX)*Mct, 1);
        end
        
        % Detect PAM signal
        function dataRX = demod(this, yd)
            dataRX = sum(bsxfun(@ge, yd, this.b.'), 2);
            dataRX = bin2gray(dataRX, 'pam', this.M).';
        end
        
        % Calculate BER in AWGN channel where the noise standard deviation 
        % is given by the function noise_std
        % noise_std is a handle function that calculates the noise std for
        % a given signal level
        function ber = ber_awgn(this, noise_std)
            ser = 0;
            for k = 1:this.M
                if k == 1
                    ser = ser + qfunc((this.b(1) - this.a(1))/noise_std(this.a(1)));
                elseif k == this.M
                    ser = ser + qfunc((this.a(k) - this.b(k-1))/noise_std(this.a(k)));
                else
                    ser = ser + qfunc((this.b(k) - this.a(k))/noise_std(this.a(k)));
                    ser = ser + qfunc((this.a(k) - this.b(k-1))/noise_std(this.a(k)));
                end
            end

            ser = ser/this.M;
            
            ber = ser/log2(this.M);
        end

        %% Level spacing (a) and decision threshold (b) optmization
        % Assumes infinite extinction ratio at first, then corrects power and
        % optmize levels again
        % The levels and thresholds calculated are after APD amplification
        function [aopt, bopt] = optimize_level_spacing_gauss_approx(this, BERtarget, rexdB, calc_noise_std, verbose)           
            % Error probability under a single tail for a given symbol
            Pe = log2(this.M)*BERtarget*(this.M/(2*(this.M-1)));

            % Initialize levels and thresholds
            aopt = zeros(this.M, 1);
            bopt = zeros(this.M-1, 1);

            rex = 10^(-abs(rexdB)/10);

            tol = Inf;
            k = 1;
            while tol(end) > this.maxtol && k < this.maxit
                apast = aopt;
                aopt(1) = aopt(end)*rex;

                for level = 1:this.M-1
                    % Find threshold
                    sig = calc_noise_std(aopt(level));

                    [dPthresh, ~, exitflag] = fzero(@(dPthresh) qfunc(abs(dPthresh)/sig) - Pe, 0);

                    if exitflag ~= 1
                        warning('level_spacing_optm: threshold optimization did not converge');
                    end

                    bopt(level) = aopt(level) + abs(dPthresh);

                    % Find next level  
                    [dPlevel, ~, exitflag] = fzero(@(dPlevel) qfunc(abs(dPlevel)/calc_noise_std(bopt(level) + abs(dPlevel))) - Pe, 0);    

                    if exitflag ~= 1
                        warning('level_spacing_optm: level optimization did not converge');     
                    end

                    aopt(level+1) = bopt(level) + abs(dPlevel);
                end

                tol(k) = sqrt(sum(abs(aopt-apast).^2));
                k = k + 1;       
            end
            
            this.b = bopt/aopt(end);
            this.a = aopt/aopt(end);

            if nargin == 5 && verbose
                figure, hold on
                plot(log(tol))
                plot([1 k], log(this.maxtol*[1 1]), 'r')
                xlabel('Iteration')
                ylabel('log(Tolerance)')
                legend('Tolerance', 'Required for Convergence')
                title('Level optimization convergece')
            end 
        end
    end
end