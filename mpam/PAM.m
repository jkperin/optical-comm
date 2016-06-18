classdef PAM
    %% Class PAM
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
        optimize_level_spacing % logical variable
    end
    
    properties (GetAccess=private)
        % Used in level spacing optimization
        maxtol = 1e-6; % maximum tolerance for convergence
        maxBERerror = 1e-3; % maximum acceptable BER error in level spacing optimization
        maxit = 20; % maximum number of iteratios
    end
       
    methods
        function obj = PAM(M, Rb, level_spacing, pshape)
            %% Class constructor
            % Inputs
            % - M = constellation size
            % - Rb = bit rate
            % - level_spacing = 'equally-spaced' or 'optimized'
            % - pshape = function handle of pulse shape 
            
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
        
        function PAMTable(self)
            %% Generate table summarizing class values
            disp('PAM class parameters summary:')
            rows = {'PAM order'; 'Symbol rate', 'Level spacing'};
            Variables = {'M'; 'Rs'; 'level_spacing'};
            Values = {self.M; self.Rb/(log2(sim.M)*1e9); self.level_spacing};
            Units = {''; 'Gbaud'; ''};

            PINtable = table(Variables, Values, Units, 'RowNames', rows)
        end
        
        %% Get methods
        function Rs = get.Rs(this)
            %% Symbol-rate assuming rectangular pulse
            Rs = this.Rb/log2(this.M);
        end
        
        function optimize_level_spacing = get.optimize_level_spacing(this)
            %% True if level_spacing == 'optimized'
            optimize_level_spacing = strcmp(this.level_spacing, 'optimized');
        end
        
        %% Auxiliary functions       
        function this = set_levels(this, levels, thresholds)
            %% Set levels to desired values
            % Levels and decision thresholds are normalized that last level is unit
            this.a = levels;
            this.b = thresholds;
        end
        
        function this = norm_levels(this)
            %% Normalize levels so that last level is unit
            this.b = this.b/this.a(end);
            this.a = this.a/this.a(end);
        end
                    
        function this = adjust_levels(this, Ptx, rexdB)
            %% Adjust levels to desired transmitted power and extinction ratio
            % Inputs:
            % - Ptx = transmitted power (W)
            % - rexdB = extinction ratio (dB). Defined as Pmin/Pmax
            % Outputs:
            % - Plevels, Pthresh = result levels and decision thresholds,
            % respectively.
            
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
        
        function [xt, xd] = mod(this, dataTX, Mct)
            %% Generate PAM signal
            % Currently, assumes that pulse shape is rectangular or that
            % there is not time overlapping between pulses
            % Inputs:
            % - dataTX = transmitted symbols from 0 to M-1
            % - Mct = oversampling rate of continuous time
            % Outputs:
            % - xt = signal in "continuous time"
            % - xd = symbols at symbol rate
            xd = this.a(gray2bin(dataTX, 'pam', this.M) + 1);
            xt = reshape(kron(xd, this.pshape(0:Mct-1)).', length(dataTX)*Mct, 1);
        end
        
        function dataRX = demod(this, yd)
            %% Demodulate PAM signal
            % Input:
            % - y = PAM signal at symbol rate
            % Output:
            % - dataRX = detected symbols
            dataRX = sum(bsxfun(@ge, yd(:), this.b.'), 2);
            dataRX = bin2gray(dataRX, 'pam', this.M).';
        end
        
        function [ber, berk] = ber_awgn(this, noise_std)
            %% Calculate BER in AWGN channel where the noise standard deviation is given by the function noise_std
            % Input:
            % - noise_std = handle function that calculates the noise std for
            % a given signal level
            % Outputs:
            % - ber = total ber
            % - berk = ber of the kth level. This actually corresponds to 
            % bertail_levels = p(error | given symbol)p(symbol)
            ser = zeros(1, this.M);
            for k = 1:this.M
                if k == 1
                    ser(k) = ser(k) + qfunc((this.b(1) - this.a(1))/noise_std(this.a(1)));
                elseif k == this.M
                    ser(k) = ser(k) + qfunc((this.a(k) - this.b(k-1))/noise_std(this.a(k)));
                else
                    ser(k) = ser(k) + qfunc((this.b(k) - this.a(k))/noise_std(this.a(k)));
                    ser(k) = ser(k) + qfunc((this.a(k) - this.b(k-1))/noise_std(this.a(k)));
                end
            end

            ser = ser/this.M;
            
            berk = ser/log2(this.M); % this corresponds to p(error | given symbol)p(symbol)
            ber = sum(berk);
        end
        
        function PreqdBm = required_power(self, N0, BERtarget)
            %% Required received power to achieve target BER for a thermal noise limited receiver with N0
            % Responsivity is assumed to be 1
            R = 1;
            Preq = (self.M-1)*sqrt(self.Rb*N0/(2*R^2*log2(self.M)))*qfuncinv(self.M*BERtarget*log2(self.M)/(2*(self.M-1)));
            PreqdBm = 10*log10(Preq/1e-3);      
            
            BER = 1/log2(self.M)*2*(self.M-1)/self.M*qfunc(sqrt(log2(self.M)/((self.M-1)^2)*Preq^2/(self.Rb*N0/2)))
            
        end
        
        function this = optimize_level_spacing_gauss_approx(this, BERtarget, rexdB, noise_std, verbose)
            %% Level spacing (a) and decision threshold (b) optmization
            % Assumes infinite extinction ratio at first, then corrects power and
            % optmize levels again
            % The levels and thresholds calculated are at the receiver
            % (i.e., after any amplification)
            % Error probability under a single tail for a given symbol
            % Inputs:
            % - BERtarget = target BER
            % - rexdB = extinction ratio (dB). Defined as Pmin/Pmax
            % - noise_std = handle function that calculates the noise std for
            % a given signal level
            % - verbose = whether to plot algorithm convergence curve
            % Outputs:
            % - aopt, bopt = optimized levels and decision thresholds, 
            % respectively.
            
            % Error probability
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
                    sig = noise_std(aopt(level));

                    [dPthresh, ~, exitflag] = fzero(@(dPthresh) qfunc(abs(dPthresh)/sig) - Pe, 0);

                    if exitflag ~= 1
                        warning('level_spacing_optm: threshold optimization did not converge');
                    end

                    bopt(level) = aopt(level) + abs(dPthresh);

                    % Find next level  
                    [dPlevel, ~, exitflag] = fzero(@(dPlevel) qfunc(abs(dPlevel)/noise_std(bopt(level) + abs(dPlevel))) - Pe, 0);    

                    if exitflag ~= 1
                        warning('level_spacing_optm: level optimization did not converge');     
                    end

                    aopt(level+1) = bopt(level) + abs(dPlevel);
                end

                tol(k) = sqrt(sum(abs(aopt-apast).^2));
                k = k + 1;       
            end
            
            this.b = bopt;
            this.a = aopt;

            BERerror = abs(this.ber_awgn(noise_std) - BERtarget)/BERtarget;
            if  BERerror > this.maxBERerror
                warning('PAM>optimize_level_spacing_gauss_approx: BER error %g greater than maximum acceptable error\n', BERerror);
            end                

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