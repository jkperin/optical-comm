classdef ofdm < handle
    properties
        prefix = 'DSB' % type of OFDM: {'DSB' (default): double sideband, 'SSB': single sideband, 'DC': DC-biased, 'ACO': asymmetrically clipped}
        Nc % # of subcarriers
        Nu % # of used subcarriers. Nu is only used to determine the oversampling ratio. In ACO-OFDM, Nu also includes the subcarriers used to trasmit 0. 
        CS % nominal constellation size
        Rb % bit rate
        Npos % positive cyclic prefix length after oversampling
        Nneg % negative cyclic prefix length after oversampling
        Pn   % power allocation
        CSn  % bit loading
        powerAllocationType = 'Levin-Campello-MA' % power allocation and bit loading type {'preemphasis': channel inversion and uniform bit loading, 'Levin-Campello-MA' (default): optimized power allocation and bit loading}
        HermitianSymmetry = true; % whether hermitian symmetry is enforced. % NOT YET USED
%         ACO = false; % whether OFDM should operate as ACO-OFDM
    end

    properties (Dependent)
        BW % one-sided bandwidth of OFDM channel
        fs % sampling rate
        Rs % symbol rate
        fc % subcarrier frequencies
        bn % # of bits per subcarrier
        ros % oversampling ratio = Nc/Nu
        K % signal attenuation factor due to clipping (K approx 1 for DC-OFDM, and K = 0.5 for ACO-OFDM)
    end
    
    properties(Hidden)
        dataTX
        dataRX
    end
    
    properties(Dependent, Hidden)
        Ncp % cyclic prefix length after oversampling
        B       % Total number of bits per symbol
        Pqam    % Subcarriers power for normalized constellations
    end
    
    properties(Constant, GetAccess=private)
        maxIterations = 100; % maximum number of iterations
        maxTol = 1e-3; % maximum tolerance
        frac_incl = 0.9999;    % fraction of energy to be included within CP length
        defaultMct = 15; % Default oversampling ratio to emulate continuous time when calculating cyclic prefix
    end
    
    methods
        function obj = ofdm(Nc, Nu, CS, Rb, prefix, powerAllocationType)
            %% Constructor
            obj.Nc = Nc;
            obj.Nu = Nu;
            obj.CS = CS;
            obj.Rb = Rb;     
            if exist('powerAllocationType', 'var')
                obj.powerAllocationType = powerAllocationType;
            else
                obj.powerAllocationType = 'Levin-Campello-MA';
            end
            
            if exist('prefix', 'var')
                obj.prefix = prefix;
                if strcmpi(obj.prefix, 'ACO') % set ACO-OFDM configuration
                    obj.aco_ofdm_config(); 
                end
            else
                obj.HermitianSymmetry = 'DSB';
            end
        end
        
        function ofdmTable = summary(self)
            %% Generate table summarizing class values
            disp('OFDM class parameters summary:')
            
            rows = {'Prefix'; 'Number of subcarriers'; 'Used subcarriers'; 'Nominal constellation size';...
                'Bit rate'; 'Cyclic prefix positive length'; 'Cyclic prefix negative length';...
                'Power allocation type'; 'Sampling rate'; 'Symbol rate'; 'Hermitian symmetry?'};

            Variables = {'prefix'; 'Nc'; 'Nu'; 'CS';...
                'Rb'; 'Npos'; 'Nneg';...
                'powerAllocationType'; 'fs'; 'Rs'; 'HermitianSymmetry'};
            Values = {self.prefix; self.Nc; self.Nu; self.CS;...
                self.Rb/1e9; self.Npos; self.Nneg;...
                self.powerAllocationType; self.fs/1e9; self.Rs/1e9; self.HermitianSymmetry};
            Units = {''; ''; ''; ''; 'Gbit/s'; ''; ''; ''; 'GHz'; 'Gbaud'; ''};

            ofdmTable = table(Variables, Values, Units, 'RowNames', rows);
        end
        
        function varargout = copy(self)
            %% Deep copy of OFDM
            for k = 1:nargout
                varargout{k} = ofdm(self.Nc, self.Nu, self.CS, self.Rb, self.prefix);
                varargout{k}.Npos = self.Npos;
                varargout{k}.Nneg = self.Nneg;
                varargout{k}.Pn = self.Pn;
                varargout{k}.CSn = self.CSn;
                varargout{k}.HermitianSymmetry = self.HermitianSymmetry;
            end
        end
        
        function s = var(self, P)
            %% OFDM signal variance under Gaussian approximation 
            if strcmpi(self.prefix, 'SSB')
                p = 1;
            else 
                p = 2;
            end
            
            if exist('P', 'var')
                s = p*sum(P);
            else
                s = p*sum(self.Pn);
            end
        end
        
        function [dc, Pmean, Pex] = dc_bias(self, Pn, rclip, Gdac, rexdB)
            %% Calculate DC bias required
            % Inputs:
            % - Pn: power allocation
            % - rclip: clipping ratio
            % - Gdac(optional, default = [1,...1]) DAC frequency response
            % at the used subcarriers frequency
            % - rexdB(optional, default = Inf) extinction ratio in dB
            % Outputs:
            % - dc: dc bias added
            % - Pmean: signal mean
            % - Pex: dc bias due to finite extinction ratio
                        
            if exist('Gdac', 'var') && not(isempty(Gdac))
                sig = sqrt(2*sum(Pn.*abs(Gdac).^2));
            else
                sig = sqrt(2*sum(Pn));
            end
            
            % DC bias and average value of each OFDM type
            if strcmpi(self.prefix, 'ACO') % ACO-OFDM
                dc = 0; 
                Pmean = sig/sqrt(2*pi);
                Pmax = Pmean + rclip*sig;
            else % DC-OFDM
                Pmean = 0;
                dc = rclip*sig;
                Pmax = 2*rclip*sig;
            end
            
            % Add additional dc bias due to non-ideal extinction ratio
            Pex = 0;
            if exist('rexdB', 'var')
                rex = 10^(-abs(rexdB)/10);
                Pex = rex*Pmax;
                dc = dc + Pex;
            end
            
            Pmean = Pmean + dc;
        end
        
        function [Padj, dc] = adjust_power_allocation(self, Pn, Ptx, rclip, rexdB)
            %% Iteratively adjust power in subcarriers so that desired average power and extinction ration are achieved
            % Inputs: 
            % - Pn: starting power allocation, which should already include
            % the attenuation due to the DAC frequency response
            % - Ptx: target average power
            % - rclip: clipping ration
            % - rexdB(optional, default = Inf): extinction ratio in dB
            
            tol = Inf;
            n = 0;
            Padj = 1;
            while tol > self.maxTol && n < self.maxIterations
                [dc, Pmean] = dc_bias(self, Padj*Pn, rclip, [], rexdB);
                Padj = Padj*(Ptx/Pmean)^2;
                tol = abs(Pmean - Ptx)/Ptx;
                n = n + 1;
            end
            
            if n == self.maxIterations
                warning('ofdm/adjust_power_allocation: max number of iterations was reached.')
            end
        end

        function power_allocation(this, Gch, varNoise, BERtarget, verbose)
            %% Power allocation and bit loading
            if strcmpi(this.prefix, 'ACO') % even subcarriers are not used in ACO-OFDM
                Gch(2:2:end) = 0; 
            end
            
            switch lower(this.powerAllocationType)
                case 'preemphasis' % channel inversion with uniform bit loading
                    [this.Pn, this.CSn] = preemphasis(this, Gch, varNoise, BERtarget);   
                case 'levin-campello-ma' % Optimal power allocation and bit loading
                    [this.Pn, this.CSn] = palloc(this, Gch, varNoise, BERtarget);
                case 'none'
                    this.Pn = ones(size(this.fc));
                    this.CSn = this.CS*ones(size(this.fc));
                otherwise 
                    error('OFDM/power_allocation: invalid option') 
            end
            
            if exist('verbose', 'var') && verbose
                figure(402), clf
                subplot(211), hold on, box on
                stem(this.fc/1e9, this.Pn/this.Pn(1), 'b', 'fill')
                plot(this.fc/1e9, abs(Gch).^2, 'r', 'LineWidth', 2)
                xlabel('Subcarrier frequency (GHz)', 'FontSize', 12)
                ylabel('Power (normalized)', 'FontSize', 12)
                title('Power allocation and bit loading')
                set(gca, 'FontSize', 12)
                grid on

                subplot(212), hold on, box on
                stem(this.fc/1e9, this.CSn, 'b', 'fill')
                xlabel('Subcarrier frequency (GHz)', 'FontSize', 12)
                ylabel('Constellation size', 'FontSize', 12)
                set(gca, 'FontSize', 12)
                ytick = sort(unique(this.CSn(this.CSn >= 4)), 'ascend');
                set(gca, 'ytick', ytick(max(end-3, 1):end))
                grid on                
                drawnow
            end
        end
        
        function [xncp, symbsTXm] = signal(this, Nsymb)
            %% Generate OFDM signal with Nsymb symbols
            % xncp : 1 x N signal with cyclic prefix
            % dataTXm : Nu/2 x Nsymb matrix containing modulated data. Each
            % column is a symbol
            % ofdm.dataTX is changed as a result of this function
            this.dataTX = cell2mat(arrayfun(@(M) randi([0 max(M-1, 0)], [1 Nsymb]), this.CSn.', 'UniformOutput', false)); % data to be encoded. Row k contains Nsymb random integers from 0 to CSn(k)
            symbsTXm = this.mod(this.dataTX); % encoded QAM symbols to be modulated onto subcarriers
            % Adjust power according to power allocation defined by Pn
            dataTXm = bsxfun(@times, symbsTXm, sqrt(this.Pn./this.Pqam).'); % Row k: sqrt(this.Pn(k)/Pqam(k))*symbsTXm(k,:)               
            
            if strcmpi(this.prefix, 'SSB') % Single sideband OFDM
                negativeSideband = zeros(size(dataTXm));
            else                
                negativeSideband = flipud(conj(dataTXm));
            end

            % zero-pad and ensure Hermitian symmetry
            % -f(N) ... -f(1) f(0) f(1) ... f(N-1) => (0 ... 0, X(-n)*, 0, X(n), 0 ... 0)
            Xn = ifftshift([zeros((this.Nc-this.Nu)/2, Nsymb); ...     % zero-pad neg. frequencies
                          negativeSideband;...                % data*(-n) or 0 (Nu) 
                          zeros(1, Nsymb); ...                % 0 at f == 0 (1)
                          dataTXm; ...                         % data(n)  (Nu)
                          zeros((this.Nc-this.Nu)/2 - 1, Nsymb)], 1);   % zero-pad pos. frequencies

            % Perform ifft (Nc-IDFT) columnwise
            xn = this.Nc*ifft(Xn, this.Nc, 1); 

            % Insert cyclic prefix           
            xncp = [xn(end-this.Ncp+1:end, :); xn]; % insert cyclic prefix

            % Parallel to serial
            xncp = reshape(xncp, 1, (this.Nc + this.Ncp)*Nsymb); % time-domain ofdm signal w/ cyclic prefix
            
            if strcmpi(this.prefix, 'ACO') % hard clip at zero
                xncp(xncp < 0) = 0;
            end 
        end
        
        function [Xn, AGCn, W] = detect(this, yk, eq, verbose)
            %% Detect OFDM signal
            % Inputs:
            % - yk: time-domain signal sampled at the symbol rate
            % - eq: struct containing equalizer parameters
            %    - mu: adaptation rate
            %    - Ntrain: training sequence length
            % Output:
            % - Xn: detected OFDM symbols
            
            yk = reshape(yk, this.Nc + this.Ncp, []);      % signal + noise

            % Remove cyclic prefix
            yn = circshift(yk(this.Npos+1:end-this.Nneg, :), -this.Nneg);  % signal + noise

            % Demodulate symbols 1/N*DFT()
            Yn = fft(yn, this.Nc, 1)/this.Nc;             % signal + noise      

            % used subcarrier amplitudes (complex conjugate subcarriers are ignored)
            Yn = Yn(1 + (1:this.Nu/2), :);            
          
            % Automatic gain control
            AGCn = sqrt(this.Pqam).'./sqrt(mean(abs(Yn).^2, 2));
            
            Yn = bsxfun(@times, Yn, AGCn);
            
            if not(isempty(eq))
                disp('OFDM: performing frequency-domain equalization...')
                % Adaptive frequency domain equalization
                if isfield(eq, 'W0') % initial guess for the weights was provided
                    if size(eq.W0, 1) < size(eq.W0, 2) % make sure this is a column vector
                        eq.W0 = eq.W0.';
                    end
                    W = eq.W0;
                else % otherwise set all weigths to 1 
                    W = ones(size(Yn, 1), 1);
                end
                Xn = zeros(size(Yn));
                e = zeros(size(Yn));
                for symb = 1:size(Yn, 2)
                    Xn(:, symb) = Yn(:, symb).*W;
                    if symb < eq.Ntrain % training sequence
                        e(:, symb) = eq.trainSeq(:, symb) - Xn(:, symb);
                        W = W + 2*eq.mu*conj(Yn(:, symb)).*e(:, symb);
                    else % decision directed
                        e(:, symb) = this.moddemod(Xn(:, symb)) - Xn(:, symb);
                        W = W + 2*eq.mu*conj(Yn(:, symb)).*e(:, symb);
                    end
                end
                % Note: 0 subcarriers in ACO-OFDM are not treated differently
                % in equalization, but their outcome is irrelevant
            else
                disp('OFDM: equalizer parameters not provided. Skipping equalization.')
                W = ones(size(Yn, 1), 1);
                Xn = Yn;
            end
            
            this.dataRX = this.demod(Xn);
            
            if exist('verbose', 'var') && verbose
                csn = find(this.CSn ~= 0, 1, 'last'); % last non-zero subcarrier
                int1 = find(this.CSn(1:ceil(this.Nu/6)) ~= 0, 1, 'last');
                int2 = ceil(this.Nu/6)+find(this.CSn(ceil(this.Nu/6)+1:ceil(this.Nu/3)) ~= 0, 1, 'last');
                idx = [1 int1 int2 csn];
                figure(403), clf
                for k = 1:4
                    subplot(2, 2, k), hold on, box on
                    plot_constellation(Xn(idx(k), :), this.dataTX(idx(k), :), this.CSn(idx(k)))
                    title(sprintf('Subcarrier: #%d', idx(k)))
                    axis square
                end
                
                if not(isempty(eq))
                    figure(404), clf, hold on, box on
                    for k = 1:4
                        subplot(2, 2, k)
                        plot(abs(e(idx(k), :)).^2)
                        xlabel('Iteration')
                        ylabel('MSE')
                        title(sprintf('Adapt. MSE for subcarrier #%d', idx(k)))
                    end

                    figure(405), clf, hold on, box on
                    plot(this.fc(1:length(W))/1e9, abs(W).^2, '-o')
                    xlabel('Frequency (GHz)')
                    ylabel('|W(f)|^2')
                    title('Equalizer')
                end
                drawnow
            end            
        end
        
        function data = demod(self, x)
            %% Recover data from symbols x
            % x must be a matrix in the form [N x Nsymb], where N denotes the number of sucarriers and Nsymb denotes the number of symbols 
            data = bsxfun(@(x, M) qamdecod(x, M), x.', self.CSn);
            data = data.';
            
            function X = qamdecod(x, M)
                % QAM decoding. Same as Matlab's qamdemod, but accounts for
                % when M = 0
                if M == 0
                    X = zeros(size(x));
                else
                    X = qamdemod(x, M, 0, 'gray');
                end
            end
        end
        
        function Xn = mod(self, data)
            %% Generate symbols from data.
            % data must be a matrix in the form [N x Nsymb], where N denotes the number of sucarriers and Nsymb denotes the number of symbols 
            Xn = bsxfun(@(data, M) qamenc(data, M), data.', self.CSn);
            Xn = Xn.';
            
            function X = qamenc(data, M)
                % QAM encoding. Same as Matlab's qammod, but accounts for
                % when M = 0
                if M == 0
                    X = zeros(size(data));
                else
                    X = qammod(data, M, 0, 'gray');
                end
            end
        end
        
        function Xn = moddemod(self, x)
            %% Equivalent to mod(demod(x)). Used to speed up decision-directed equalization
            % x must be a column vector
            Xn = arrayfun(@(x, M) qamencdec(x, M), x, self.CSn.'); % Row k: qamencdec(x(k), CSn(k))
            
            function X = qamencdec(x, M)
                % QAM encoding. Same as qammod(qamdemod()), but accounts for
                % when M = 0
                if M == 0
                    X = 0;
                else
                    X = qammod(qamdemod(x, M, 0, 'gray'), M, 0, 'gray');
                end
            end
        end
        
        function [bercount, interval] = count_ber(self, Ndiscard, verbose)
            %% Measured BER
            ind = Ndiscard(1)+1:size(self.dataTX, 2)-Ndiscard(2);
            Nsymb = length(ind);
            Nerr = zeros(1, self.Nu/2);
            Nbits = 0;
            for kk = 1:size(self.dataRX, 1)
                if self.CSn(kk) ~= 0                    
                    Nerr(kk) = biterr(self.dataTX(kk, ind), self.dataRX(kk, ind), self.bn(kk));
                    Nbits = Nbits + Nsymb*self.bn(kk);
                end
            end

            % average BER
            bercount = sum(Nerr)/Nbits;

            % 95% confidence intervals for the counted BER
            [~, interval] = berconfint(Nerr, Nbits);
            interval = [min(interval(:,1)), max(interval(:,2))];  
            
            if exist('verbose', 'var') && verbose
                figure(410), clf, box on
                stem(1:self.Nu/2, Nerr, 'fill')
                xlabel('Subcarrier')
                ylabel('Number of bit errors')
                drawnow
            end
        end
               
        function ber_theory = calc_ber(self, SNRdBn)
            %% Calculate theoretical BER from SNR in dB at each subcarrier
            ber = zeros(1, length(SNRdBn));
            for k = 1:length(SNRdBn)
                if self.CSn(k) ~= 0
                    ber(k) = berqam(self.CSn(k), SNRdBn(k));
                end
            end

            ber_theory = sum(ber.*self.bn)/sum(self.bn);
        end
        
        function [ber_est, SNRn] = estimate_ber(this, Pnmeasured, Pnoise, Gch, noiseVar, verbose)
            %% Estimate BER from signal
            % Inputs:
            % - Pnmeasured: estimate of the power levels in each subcarrier.
            % Contains signal + noise
            % - Gch: channel frequency response at the subcarriers frequency
            % - noiseVar: noise variance at the subcarriers frequency
            % - verbose (optional, default = false): whether to plot
            % results
            Gch = this.K*Gch; % accounts for attenuation due to clipping
            Pnest = this.Pn.*abs(Gch).^2; % power at subcarriers at the receiver
            
            % Calculate unbiased SNR estimate
            SNRn.measured = 10*log10(Pnmeasured./Pnoise - 1); 
            SNRn.estimated = 10*log10(Pnest./noiseVar - 1); 
            
            ber_est = this.calc_ber(SNRn.estimated);
            
            if exist('verbose', 'var') && verbose
                idx = find(this.CSn ~= 0);
                fcGHz = this.fc(idx)/1e9;
                figure(301), clf
                subplot(211), hold on, box on
                stem(fcGHz, Pnmeasured(idx), 'fill')
                stem(fcGHz, Pnest(idx))
                xlabel('Frequency (GHz)')
                ylabel('Received subcarrier power')
                legend('Measured (signal + noise)', 'Estimated (signal)');
                subplot(212), hold on, box on
                stem(fcGHz, Pnoise(idx), 'fill')
                stem(fcGHz, noiseVar(idx))
                xlabel('Frequency (GHz)')
                ylabel('Noise variance')
                legend('Measured', 'Estimated')
                figure(302), clf, hold on, box on
                stem(fcGHz, SNRn.measured(idx), 'fill')
                stem(fcGHz, SNRn.estimated(idx))
                xlabel('Frequency (GHz)')
                ylabel('Unbiased SNR (dB)')    
                legend('Measured', 'Estimated')
                drawnow
            end
        end
        
        function varQ = quantization_noise_var(self, ENOB, rclip)
            %% Quantization noise variance as a function of the signal power
            if isinf(ENOB)
                varQ = @(Psig) 0;
            else
                % The general expression for quantization noise is
                % varQ = (1-Pc)/12*(DeltaX/(2^ENOB-1))^2;
                % This assumes that there's no quantization error at the
                % clipping levels.
                if strcmpi(self.prefix, 'ACO') % ACO-OFDM
                    Pc = 1/2 + qfunc(rclip); % clipping probability
                    varQ = @(Psig) (1-Pc)/12*(Psig*sqrt(2*pi)*rclip/(2^ENOB))^2; % signal goes from 0 to rclip*sig
                else % DC-OFDM
                    Pc = 2*qfunc(rclip); % clipping probability
                    varQ = @(Psig) (1-Pc)/12*(2*Psig/(2^ENOB))^2; % signal goes from 0 to 2Psig (two times the average power, since clipping in DC-OFDM is symmetric)
                end                
            end
        end
        
        function set_cyclic_prefix(self, Nn, Np)
            %% Set cyclic prefix length
            self.Nneg = Nn;
            self.Npos = Np;
        end
               
        function [Nn, Np] = cyclic_prefix(this, Hch, verbose)
        %% Calculate the cyclic prefix length after oversampling 
        % Inputs:
        % - Hch: function handle of all input responses of the channel combined
        % Hch(f, fs) where f is the frequency vector and fs is the sampling
        % rate
        % Algorithm:
        % 1. Assumes a cyclic prefix length k
        % 2. From k, it calculates the new sampling rate fs
        % 3. Calculate the number of samples at both sides (Nneg and Npos) that 
        % contains the desired fraction of energy
        % 4. From that calculates this.Ncp = Npos + Neg; (zero is not included)
        % 5. If this.Ncp == k end simulation, otherwise increment k and repeat

%             this.Ncp = 0;
            Mct = this.defaultMct;
            N = 2^14;
            
            % Time and frequency 
            fsamp = this.fs*Mct;
            [f, t] = freq_time(N, fsamp);
            t = t - 0.5*N/fsamp;
            
            % Channel frequency response
            Hch = Hch(f, fsamp);
            Hch = Hgrpdelay(Hch, f); % Calculate and remove group delay
            
            % Time domain
            ht = fftshift(real(ifft(ifftshift(Hch))));
            ht = ht/max(ht);
            
            % Initialization
            Nd = ceil(this.Nc/this.ros);
            k = 0;
            Npre = Inf; % cyclic prefix length
            n = -512:512;
            while Npre > k && k < this.maxIterations
                fs_temp = this.Rs*(this.Nc + k)/Nd;
                
                tn = n/fs_temp;
                hn = interp1(t, ht, tn); % retime

                % CP based on energy
                en_frac = cumsum(hn.^2)/sum(hn.^2);
                Nn = sum(en_frac >= (1 - this.frac_incl)/2 & n < 0);
                Np = sum(en_frac <= (1 + this.frac_incl)/2 & n > 0);

                Npre = Nn + Np;
                
                k = k + 1;
            end         

            assert(k ~= this.maxIterations, 'ofdm/cyclic_prefix: CP calculation did not converge');
            
            if isempty(this.Nneg)
                this.Nneg = Nn;
                this.Npos = Np;
            else
                if this.Nneg ~= Nn || this.Npos ~= Np
                    warning('ofdm/cyclic_prefix: cyclic prefix length is not equal to the desired length.')
                    fprintf('Nneq(set) = %d, Nneq(required) = %d\n', this.Nneg, Nn)
                    fprintf('Npos(set) = %d, Npos(required) = %d\n', this.Npos, Np)
                end
            end
            
            if exist('verbose', 'var') && verbose
                figure(401), grid on, hold on, box on
                plot(t, ht)
                stem(tn, hn, 'fill')
                plot(-this.Nneg/this.fs*[1 1], [-1 1], 'k')
                plot(this.Npos/this.fs*[1 1], [-1 1], 'k')
                plot(-Nn/this.fs*[1 1], [-1 1], ':k')
                plot(Np/this.fs*[1 1], [-1 1], ':k')
                text(-(this.Nneg)/this.fs, 0.7, sprintf('N_{neg} = %d', this.Nneg))
                text((this.Npos+0.1)/this.fs, 0.7, sprintf('N_{pos} = %d', this.Npos))
                axis([1.5*this.Ncp*[-1 1]*1/this.fs -0.5 1])
                xlabel('t (s)')
                ylabel('p(t)')
                title('Impulse response of the channel')
                drawnow
            end            
        end       
    end
    
    %% Get and Set Methods
    methods 
        function fs = get.fs(this)
            %% OFDM sampling rate
            fs = this.Rs*this.ros*(this.Nc + this.Ncp)/this.Nc;
        end        
        
        function ros = get.ros(this)
            %% OFDM oversampling ratio
            ros = this.Nc/this.Nu;
        end
        
        function fc = get.fc(this)
            %% Frequency at which subcarriers are located
            fc = this.fs/this.Nc*(1:this.Nu/2);   
            % Note: this does not include last subcarriers set to 0 for
            % oversampling. fc does include, however, subcarriers set to 0
            % to accomplish ACO-OFDM.
        end
        
        function BW = get.BW(self)
            %% One-sided bandwidth of OFDM signal
            if strcmpi(self.prefix, 'ACO') % ACO-OFDM
                BW = 2*self.Rb/log2(self.CS);
            else
                BW = self.Rb/log2(self.CS);
            end
        end
                
        function Rs = get.Rs(this)
            %% Symbol rate
            p = 1;
            if strcmpi(this.prefix, 'ACO') % ACO-OFDM requires twice the symbol rate
                p = 2;              
            end
            Rs = 2*p*this.Rb/log2(this.CS);
        end
       
        function B = get.B(this)
            %% Total number of bits
            if isempty(this.CSn)
                B = this.Nu/2*log2(this.CS); 
            else
                B = sum(this.bn);
            end
        end
        
        function Ncp = get.Ncp(this)
            %% Cyclic prefix length
            if any(isempty([this.Npos this.Nneg]))
                Ncp = 0;
                warning('ofdm/Ncp: Cyclic prefix length was not assigned/calculated yet. Assuming no cyclic prefix.')
            else
                Ncp = this.Npos + this.Nneg;
            end
        end
        
        function aco_ofdm_config(self)
            %% Puts OFDM class in ACO-OFDM configure, whereby zero subcarriers are set to 0
            self.prefix = 'ACO';
            self.CSn = zeros(1, self.Nu/2);
            self.CSn(1:2:end) = self.CS; % even subcarriers are not used.
        end
            
        function Pqam = get.Pqam(this)
        	%% Calculate the average power of a CS-QAM constellation with dmin = 2 (used to normalize symbol power after qammod)
            Pqam = ones(1, this.Nu/2);
            for k = 1:this.Nu/2
                if this.CSn(k) == 0
                    continue
                else
                    Pqam(k) = mean(abs(qammod(0:this.CSn(k)-1, this.CSn(k), 0, 'gray')).^2);
                end
            end
        end  
        
        function bn = get.bn(self)
            %% Bits per subcarrier
            bn = zeros(size(self.CSn));
            bn(self.CSn ~= 0) = log2(self.CSn(self.CSn ~= 0));
        end
        
        function set.powerAllocationType(self, type)
            %% Power allocation type
            if strcmpi(type, 'preemphasis') || strcmpi(type, 'Levin-Campello-MA') || strcmpi(type, 'none')
                self.powerAllocationType = lower(type);
            else
                error('ofdm/powerAllocationType: Invalid power allocation type. It must be either preemphasis, Levin-Campello-MA, or none.')
            end
        end
        
        function Pn = get.Pn(self)
            %% Get Pn
            if isempty(self.Pn)
                Pn = ones(1, self.Nu/2);
                if strcmpi(self.prefix, 'ACO')
                    Pn(2:2:end) = 0;
                end
            else
                Pn = self.Pn;
            end
        end
        
        function CSn = get.CSn(self)
            %% Get CSn
            if isempty(self.CSn)
                CSn = self.CS*ones(1, self.Nu/2);
                if strcmpi(self.prefix, 'ACO')
                    CSn(2:2:end) = 0;
                end
            else
                CSn = self.CSn;
            end
        end
        
        function K = get.K(self)
            %% Signal attenuation due to clipping
            % After clipping, xc(t) = Kx(t) + d(t). Hence subcarriers 
            % power is scaled by K^2. x(t) and d(t) are uncorrelated. This follows
            % from Bussgang's theorem of nonlinearities in Gaussian-distributed
            % signals
            % K = 1 - Q(r+) - Q(r-), where Q() is qfunc, and r+ and r-
            % denote the positive and negative clipping ratios
            if strcmpi(self.prefix, 'ACO')
                K = 0.5;
            else
                K = 1;
            end
        end
        
        function clear(self)
            %% Clear heavy variables
            self.dataTX = [];
            self.dataRX = [];
        end
    end
end