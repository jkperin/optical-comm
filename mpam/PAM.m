classdef PAM
    %% Class PAM
    properties
        M % constellation size
        Rb % bit rate
        level_spacing % 'equally-spaced' or 'optimized'
        pulse_shape % struct containing the following fields {type: {'rectangular', 'root raised cosine', 'raised cosine'}, h: pulse shape impulse response i.e., coefficients of FIR filter,sbs: samples per symbol, and other parameters such as rolloff factor, etc}
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
        maxit = 50; % maximum number of iteratios
    end
       
    methods
        function obj = PAM(M, Rb, level_spacing, pulse_shape)
            %% Class constructor
            % Inputs
            % - M = constellation size
            % - Rb = bit rate
            % - level_spacing = 'equally-spaced' or 'optimized'
            % - puse_shape = % struct containing the following fields {type:
            % {'rectangular', 'root raised cosine', 'raised cosine'}, 
            % h: pulse shape impulse response i.e., coefficients of FIR filter,
            % Mct: oversampling ratio of pulse shapping filter, and
            % other parameters such as rolloff factor, etc}
            
            obj.M = M;
            obj.Rb = Rb;
            
            if exist('level_spacing', 'var')
                obj.level_spacing = level_spacing;
            else
                obj.level_spacing = 'equally-spaced';
            end
            
            if exist('pulse_shape', 'var')
                assert(~strcmp(class(pulse_shape), 'function_handle'), 'PAM: behaviour of PAM has changed. Now PAM expects a struct as opposed to the function handle.')
                obj.pulse_shape = pulse_shape;
                obj.pulse_shape.h = obj.norm_filter_coefficients(obj.pulse_shape.h);
            else
                obj.pulse_shape = select_pulse_shape('rect', 1);
            end
            
            obj = obj.reset_levels();
        end
        
        function PAMtable = summary(self)
            %% Generate table summarizing class values
            disp('PAM class parameters summary:')
            rows = {'PAM order'; 'Symbol rate'; 'Level spacing'; 'Pulse shape'; 'Samples per symbol in pulse shaping'};
            Variables = {'M'; 'Rs'; 'level_spacing'; 'pulse_shape.type'; 'pulse_shape.sps'};
            Values = {self.M; self.Rs; self.level_spacing; self.pulse_shape.type; self.pulse_shape.sps};
            Units = {''; 'Gbaud'; ''; ''; ''};

            PAMtable = table(Variables, Values, Units, 'RowNames', rows);
        end
    end
    
    methods 
        %% Get and set methods
        function Rs = get.Rs(self)
            %% Symbol-rate assuming rectangular pulse
            Rs = self.Rb/log2(self.M);
        end
        
        function optimize_level_spacing = get.optimize_level_spacing(self)
            %% True if level_spacing == 'optimized'
            optimize_level_spacing = strcmp(self.level_spacing, 'optimized');
        end
        
        function H = Hpshape(self, f)
            %% Frequency response of PAM pulse shape
            fs = self.Rs*self.pulse_shape.sps;
            delay = grpdelay(self.pulse_shape.h, 1, 1);
            H = freqz(self.pulse_shape.h/abs(sum(self.pulse_shape.h)), 1, f, fs)...
                .*exp(1j*2*pi*f/fs.*delay); % remove group delay
        end        
    end
       
    methods
        %% Levels and decision thresholds     
        function self = set_levels(self, levels, thresholds)
            %% Set levels to desired values
            % Levels and decision thresholds are normalized that last level is unit
            assert(length(levels) == self.M, 'mpam/set_levels: invalid number of levels');
            assert(length(thresholds) == self.M-1, 'mpam/set_levels: invalid number of decision thresholds');
            if size(levels, 1) < size(levels, 2) % ensure levels are M x 1 vector
                levels = levels.'; 
            end
            if size(thresholds, 1) < size(thresholds, 2)  % ensure thresholds are M x 1 vector
                thresholds = thresholds.';
            end
                
            self.a = levels;
            self.b = thresholds;
        end
        
        function self = norm_levels(self)
            %% Normalize levels so that last level is 1
            self.b = self.b/self.a(end);
            self.a = self.a/self.a(end);
        end
        
        function self = unbias(self)
            %% Remove DC bias from levels and normalize to have excusion from -1 to 1
            self.b = self.b - mean(self.a);
            self.a = self.a - mean(self.a);
            self = self.norm_levels;
        end
            
        function self = reset_levels(self)
            %% Reset levels and decision thresholds to original configuration
            self.a = ((0:2:2*(self.M-1))/(2*(self.M-1))).';
            self.b = ((1:2:(2*(self.M-1)-1))/(2*(self.M-1))).';
        end
                    
        function self = adjust_levels(self, Ptx, rexdB)
            %% Adjust levels to desired transmitted power and extinction ratio
            % Inputs:
            % - Ptx = transmitted power (W)
            % - rexdB = extinction ratio (dB). Defined as Pmin/Pmax
            % Outputs:
            % - Plevels, Pthresh = result levels and decision thresholds,
            % respectively.
            
            rex = 10^(-abs(rexdB)/10); % extinction ratio. Defined as Pmin/Pmax
            switch self.level_spacing
                case 'equally-spaced'
                    % Restart levels
                    self.a = ((0:2:2*(self.M-1))/(2*(self.M-1))).';
                    self.b = ((1:2:(2*(self.M-1)-1))/(2*(self.M-1))).';
                    
                    amean = mean(self.a); % mean value
                    
                    Pmin = 2*Ptx*rex/(1 + rex); % power of the lowest level 
                    Plevels = self.a*(Ptx/amean)*((1-rex)/(1+rex)) + Pmin; % levels at the transmitter
                    Pthresh = self.b*(Ptx/amean)*((1-rex)/(1+rex)) + Pmin; % decision thresholds at the transmitter
                case 'optimized'
                    amean = mean(self.a); % mean value
                    
                    % Extinction ratio was already enforced in the
                    % optimization process, so just scale to desired power
                    Plevels = self.a*Ptx/amean; % levels at the transmitter
                    Pthresh = self.b*Ptx/amean; % decision thresholds at the transmitter
                otherwise
                    error('pam class: Invalid level spacing option')
            end
            
            self.a = Plevels;
            self.b = Pthresh;            
        end   
        
        function bset = place_thresholds(self)
            %% Given leves mpam.a, find best position of thresholds by simply scaling and/or offseting them
            %% This assumes that noise is signal independent
            tideal = diff(self.a);
            tideal = self.a(1:end-1) + tideal/2;
            cost_function = @(t) norm(t - tideal);
            
            [Vset, fval, exitflag] = fminsearch(@(V) cost_function(V(1)*(self.b - mean(self.b)) + V(2)), [1 0]);
    
            if exitflag ~= 1
               warning('PAM/place_thresholds: optimization did not converge')
            end   
            
            bset = Vset(1)*(self.b - mean(self.b)) + Vset(2);
        end
        
        function self = mzm_predistortion(self, Vswing, Vbias, verbose)
            %% Predistort levels to compensate for MZM nonlinear response in IM-DD
            predist = @(p) 2/pi*asin(sqrt(p));
            dist = @(v) sin(pi/2*v)^2;
            
            Vmin = Vbias - Vswing/2;
            Vmax = Vbias + Vswing/2;
            
            Pmax = dist(Vmax);
            Pmin = dist(Vmin);
            DP = (Pmax-Pmin)/(self.M-1);
            Pk = Pmin:DP:Pmax;
            
            assert(all(Pk >= 0), 'mpam/mzm_predist: Pswing too high. Pswing must be such that all levels are non-negative')
            
            % Predistortion
            Vk = predist(Pk);
            
            % Set
            self = self.set_levels(Vk, self.b);
            % Note: Decision thresholds are not predistorted, since the predistortion is only
            % used for generatign the levels at the transmitter.
            
            if exist('verbose', 'var') && verbose
                figure(233), clf, hold on, box on
                t = linspace(0, 1);
                plot(t, sin(pi/2*t).^2, 'k');
                plot((self.a*[1 1]).', [zeros(1, self.M); sin(pi/2*self.a.').^2], 'k');
                plot([zeros(1, self.M); self.a.'], ([1; 1]*Pk), 'k')
                xlabel('Driving signal')
                ylabel('Resulting power levels')
                axis([0 1 0 1])
            end        
        end
    end
    
    methods
        %% Modulation and demodulation
        function [xt, xd] = signal(self, dataTX)
            %% Generate PAM signal
            % Note: group delay due to pulse shaping filtering is not
            % removed
            % Input:
            % - dataTX = transmitted symbols from 0 to M-1
            % Outputs:
            % - xt = pulse shaped PAM signal with oversampling ratio = mpam.pulse_shape.sps (1 x length(dataTX)*mpam.pulse_shape.sps)
            % - xd = symbols at symbol rate (1 x length(dataTX))
            
            % Normalize filter taps to preserve levels amplitude           
            self.pulse_shape.h = self.norm_filter_coefficients(self.pulse_shape.h);
            % Generate data, upsample, filter, and remove group delay
            xd = self.mod(dataTX); % 1 x length(dataTX)
            ximp = upsample(xd, self.pulse_shape.sps);
            xt = filter(self.pulse_shape.h, 1, ximp); 
        end
        
        function xd = mod(self, dataTX)
            %% Generate PAM symbols
            % Input:
            % - dataTX = transmitted symbols from 0 to M-1
            % Outputs:
            % - xd = symbols at symbol rate (1 x length(dataTX))            
            xd = self.a(gray2bin(dataTX, 'pam', self.M) + 1).'; % 1 x length(dataTX)
        end
        
        function dataRX = demod(self, yd)
            %% Demodulate PAM signal
            % Input:
            % - y = PAM signal at symbol rate
            % Output:
            % - dataRX = detected symbols
            dataRX = sum(bsxfun(@ge, yd(:), self.b.'), 2);
            dataRX = bin2gray(dataRX, 'pam', self.M).';
        end
        
        function [dataRX, mpamOpt] = demod_sweeping_thresholds(self, yd, dataTX, validInd, verbose)
            %% Demdulate PAM signal by sweeping thresholds until BER is minimized
            if not(exist('validInd', 'var'))
                validInd = 1:length(dataTX);
            end

            mpamOpt = self;
            if exist('verbose', 'var') && verbose
                figure(999), clf, hold on, box on
            end
            for k = 1:self.M-1
                [topt, fval, exitflag] = fminbnd(@(t) calc_log10_ber(t, k), mpamOpt.a(k), mpamOpt.a(k+1));
                
                if exitflag ~= 1
                    warning('PAM/demod_swiping_thresholds: Threshold swiping did not cpnverge. Threshold swiping of threshold %d exited with exitflag %d', k, exitflag);
                end
                
                ts = linspace(mpamOpt.a(k), mpamOpt.a(k+1));
                lber = zeros(size(ts));
                for kk = 1:length(ts)
                    lber(kk) = calc_log10_ber(ts(kk), k);
                end
                
                if exist('verbose', 'var') && verbose
                    figure(999)
                    plot(ts, lber)
                    plot(topt, fval, 'ok')
                    drawnow
                end
                
                mpamOpt.b(k) = topt;
            end
               
            dataRX = mpamOpt.demod(yd);
            
            function log10_ber = calc_log10_ber(threshold, threshold_idx)
                mpamOpt.b(threshold_idx) = threshold;
                
                data = mpamOpt.demod(yd);
                
                [~, ber] = biterr(dataTX(validInd), data(validInd));
                
                log10_ber = log10(ber);
            end
        end        
    end
    
    methods
        function [ber, berk] = berAWGN(self, noise_std)
            %% Calculate BER in AWGN channel where the noise standard deviation is given by the function noise_std
            % Input:
            % - noise_std = handle function that calculates the noise std for
            % a given signal level
            % Outputs:
            % - ber = total ber
            % - berk = ber of the kth level. self actually corresponds to 
            % bertail_levels = p(error | given symbol)p(symbol)
            ser = zeros(1, self.M);
            for k = 1:self.M
                if k == 1
                    ser(k) = ser(k) + qfunc((self.b(1) - self.a(1))/noise_std(self.a(1)));
                elseif k == self.M
                    ser(k) = ser(k) + qfunc((self.a(k) - self.b(k-1))/noise_std(self.a(k)));
                else
                    ser(k) = ser(k) + qfunc((self.b(k) - self.a(k))/noise_std(self.a(k)));
                    ser(k) = ser(k) + qfunc((self.a(k) - self.b(k-1))/noise_std(self.a(k)));
                end
            end

            ser = ser/self.M;
            
            berk = ser/log2(self.M); % self corresponds to p(error | given symbol)p(symbol)
            ber = sum(berk);
        end
        
        function [PreqdBm, BER] = required_power(self, N0, BERtarget)
            %% Required received power to achieve target BER for a thermal noise limited receiver with N0
            % Responsivity is assumed to be 1
            R = 1;
            Preq = (self.M-1)*sqrt(self.Rb*N0/(2*R^2*log2(self.M)))*qfuncinv(self.M*BERtarget*log2(self.M)/(2*(self.M-1)));
            PreqdBm = 10*log10(Preq/1e-3);      
            
            % Show resulting BER just for sanity check
            BER = 1/log2(self.M)*2*(self.M-1)/self.M*qfunc(sqrt(log2(self.M)/((self.M-1)^2)*Preq^2/(self.Rb*N0/2))) 
        end
        
        function self = optimize_level_spacing_gauss_approx(self, BERtarget, rexdB, noise_std, verbose)
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
            Pe = log2(self.M)*BERtarget*(self.M/(2*(self.M-1)));

            % Initialize levels and thresholds
            aopt = zeros(self.M, 1);
            bopt = zeros(self.M-1, 1);

            rex = 10^(-abs(rexdB)/10);

            tol = Inf;
            k = 1;
            while tol(end) > self.maxtol && k < self.maxit
                apast = aopt;
                aopt(1) = aopt(end)*rex;

                for level = 1:self.M-1
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
            
            self.b = bopt;
            self.a = aopt;

            BERerror = abs(self.ber_awgn(noise_std) - BERtarget)/BERtarget;
            if  BERerror > self.maxBERerror
                warning('PAM>optimize_level_spacing_gauss_approx: BER error %g greater than maximum acceptable error\n', BERerror);
            end                

            if nargin == 5 && verbose
                figure(645), hold on
                plot(log(tol))
                plot([1 k], log(self.maxtol*[1 1]), 'r')
                xlabel('Iteration')
                ylabel('log(Tolerance)')
                legend('Tolerance', 'Required for Convergence')
                title('Level optimization convergece')
                drawnow
            end 
        end           
    end
    
    %% Auxiliary functions
    methods
        function h = norm_filter_coefficients(~, h)
            %% Normalize coefficients of FIR filter h so that impulse response of h at t = 0 is 1
            n = length(h);
            if mod(n, 2) == 0 % even
                h = 2*h/(h(n/2) + h(n/2+1));
            else
                h = h/h((n+1)/2);
            end
        end
    end
            
    %% Validation methods
    methods
        function validate_pulse_shape(self)
            Mct = self.pulse_shape.sps;
            dataTX = randi([0 self.M-1], [1 1024]);
            xt = self.signal(dataTX);
            f = freq_time(length(xt), 1);
            
            figure
            eyediagram(xt, 2*Mct)
            title('Eye diagram')
            figure, box on
            plot(w/(2*pi), abs(self.Hpshape(f)).^2)
            xlabel('Normalized frequency (Hz)')
            ylabel('Amplitude')
            title('Frequency response of pulse shaping filter')
        end
    end
end