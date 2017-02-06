classdef ofdm < handle
    properties
        Nc % number of subcarriers
        Nu % number of used subcarriers 
        CS % nominal constellation size
        Rb % bit rate
        Npos % positive cyclic prefix length after oversampling
        Nneg % negative cyclic prefix length after oversampling
        Pn   % power allocation
        CSn  % bit loading
        powerAllocationType % power allocation and bit loading type {'preemphasis': channel inversion and uniform bit loading, 'palloc' (default): optimized power allocation and bit loading}
    end

    properties (Dependent)
        fs % sampling rate
        Rs % symbol rate
        Ms % OFDM oversampling ratio
        fc % subcarrier frequencies
        bn % # of bits per subcarrier
    end
    
    properties(Hidden)
        varQtx
        varQrx
        dataTX
        dataRX
    end
    
    properties(Dependent, Hidden)
        Npre_os
        B       % Total number of bits per symbol
        Pqam    % Average power of generated OFDM symbol with qammod
    end
    
    properties(Constant, GetAccess=private)
        max_iterations = 100; % maximum number of iterations
        frac_incl = 0.9999;    % fraction of energy to be included within CP length
        defaultMct = 15; % Default oversampling ratio to emulate continuous time when calculating cyclic prefix
    end
    
    methods
        function obj = ofdm(Nc, Nu, CS, Rb, powerAllocationType)
            %% Constructor
            obj.Nc = Nc;
            obj.Nu = Nu;
            obj.CS = CS;
            obj.Rb = Rb;     
            if exist('powerAllocationType', 'var')
                obj.powerAllocationType = powerAllocationType;
            else
                obj.powerAllocationType = 'palloc';
            end
        end
        
        function ofdmTable = summary(self)
            %% Generate table summarizing class values
            disp('OFDM class parameters summary:')
            
            rows = {'Number of subcarriers'; 'Used subcarriers'; 'Nominal constellation size';...
                'Bit rate'; 'Cyclic prefix positive length'; 'Cyclic prefix negative length';...
                'Power allocation type'; 'Sampling rate'; 'Symbol rate'};

            Variables = {'Nc'; 'Nu'; 'CS';...
                'Rb'; 'Npos'; 'Nneg';...
                'powerAllocationType', 'fs', 'Rs'};
            Values = {self.Nc; self.Nu; self.CS;...
                self.Rb/1e9; self.Npos; self.Nneg;...
                self.powerAllocationType; self.fs/1e9; self.Rs/1e9};
            Units = {''; ''; ''; 'Gbit/s'; ''; ''; ''; 'GHz', 'Gbaud'};

            ofdmTable = table(Variables, Values, Units, 'RowNames', rows);
        end
        
        
        function s = var(self, P)
            %% OFDM signal variance under Gaussian approximation 
            if exist('P', 'var')
                s = 2*sum(P);
            else
                s = 2*sum(self.Pn);
            end
        end
        
        function dc = dc_bias(~, Pn, rclip, Gdac, rexdB)
            %% Calculate DC bias required for OFDM
            if exist('Gdac', 'var') && not(isempty(Gdac))
                dc = rclip*sqrt(2*sum(Pn.*abs(Gdac).^2));
            else
                dc = rclip*sqrt(2*sum(Pn));
            end
            
            % Add additional dc bias due to non-ideal extinction ratio
            if exist('rexdB', 'var')
                rex = 10^(-abs(rexdB)/10);
                mindc = rex*2*dc;
                dc = dc + mindc;
            end
        end

        %% Power allocation and bit loading
        function power_allocation(this, Gch, varNoise, BERtarget, verbose)
            %% Preemphasis at the transmitter
            switch this.powerAllocationType
                case 'preemphasis'
                    [this.Pn, this.CSn] = preemphasis(this, Gch, varNoise, BERtarget);
                %% Optimal power allocation and bit loading   
                case 'palloc'
                    [this.Pn, this.CSn] = palloc(this, Gch, varNoise, BERtarget);
                case 'none'
                    this.Pn = ones(size(this.fc));
                    this.CSn = this.CS*ones(size(this.fc));
                otherwise 
                    error('power_allocation: invalid option') 
            end
            
            if exist('verbose', 'var') && verbose
                figure(402), clf
                subplot(211), hold on, box on
                stem(this.fc/1e9, this.Pn/this.Pn(1))
                plot(this.fc/1e9, abs(Gch).^2, 'LineWidth', 2)
                xlabel('Subcarrier frequency (GHz)', 'FontSize', 12)
                ylabel('Power (normalized)', 'FontSize', 12)
                title('Power allocation and bit loading')
                set(gca, 'FontSize', 12)
                grid on

                subplot(212), hold on, box on
                stem(this.fc/1e9, this.CSn)
                xlabel('Subcarrier frequency (GHz)', 'FontSize', 12)
                ylabel('Constellation size', 'FontSize', 12)
                set(gca, 'FontSize', 12)
                set(gca, 'ytick', [4 16 32 64])
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
            this.dataTX = zeros(this.Nu/2, Nsymb);
            dataTXm = zeros(this.Nu/2, Nsymb);
            symbsTXm = zeros(this.Nu/2, Nsymb);
            for kk = 1:this.Nu/2
                if this.CSn(kk) > 1
                    this.dataTX(kk,:) = randi([0 this.CSn(kk)-1], [1 Nsymb]);              % data to be encoded (symbols columnwise)
                    symbsTXm(kk,:) = qammod(this.dataTX(kk,:), this.CSn(kk), 0, 'gray');    % encoded QAM symbols to be modulated onto subcarriers
                    dataTXm(kk,:) = sqrt(this.Pn(kk))*symbsTXm(kk,:)/sqrt(this.Pqam(kk));  % scale constellation so that Pn is the power at nth subcarrier (E(|Xn|^2) = Pn)
                else
                    continue;
                end
            end          
            
            % zero-pad and ensure Hermitian symmetry
            % -f(N) ... -f(1) f(0) f(1) ... f(N-1) => (0 ... 0, X(-n)*, 0, X(n), 0 ... 0)
            Xn = ifftshift([zeros((this.Nc-this.Nu)/2, Nsymb); ...     % zero-pad neg. frequencies
                          flipud(conj(dataTXm));...            % data*(-n) (Nu) 
                          zeros(1, Nsymb); ...                % 0 at f == 0 (1)
                          dataTXm; ...                         % data(n)  (Nu)
                          zeros((this.Nc-this.Nu)/2 - 1, Nsymb)], 1);   % zero-pad pos. frequencies

            % Perform ifft (Nc-IDFT) columnwise
            xn = this.Nc*ifft(Xn, this.Nc, 1); 

            % Insert cyclic prefix           
            xncp = [xn(end-this.Npre_os+1:end, :); xn]; % insert cyclic prefix

            % Parallel to serial
            xncp = reshape(xncp, 1, (this.Nc + this.Npre_os)*Nsymb); % time-domain ofdm signal w/ cyclic prefix
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
            
            yk = reshape(yk, this.Nc + this.Npre_os, []);      % signal + noise

            % Remove cyclic prefix
            yn = circshift(yk(this.Npos+1:end-this.Nneg, :), -this.Nneg);  % signal + noise

            % Demodulate symbols 1/N*DFT()
            Yn = fft(yn, this.Nc, 1)/this.Nc;             % signal + noise      

            % used subcarrier amplitudes (complex conjugate subcarriers are ignored)
            Yn = Yn(1 + (1:this.Nu/2), :);            
          
            % Automatic gain control
            AGCn = sqrt(this.Pqam).'./sqrt(mean(abs(Yn).^2, 2));
            
            Yn = bsxfun(@times, Yn, AGCn);
            
            % Adaptive frequency domain equalization
            W = ones(size(Yn, 1), 1);
            Xn = zeros(size(Yn));
            e = zeros(size(Yn));
            for symb = 1:size(Yn, 2)
                Xn(:, symb) = Yn(:, symb).*W;
                if symb < eq.Ntrain % training sequence
                    e(:, symb) = eq.trainSeq(:, symb) - Xn(:, symb);
                    W = W + 2*eq.mu*conj(Yn(:, symb)).*e(:, symb);
                else % decision directed
                    e(:, symb) = this.mod(this.demod(Xn(:, symb))) - Xn(:, symb);
                    W = W + 2*eq.mu*conj(Yn(:, symb)).*e(:, symb);
                end
            end
            
            this.dataRX = this.demod(Xn);
            
            if exist('verbose', 'var') && verbose
                figure(403), clf
                subplot(221), hold on, box on
                plot_constellation(Xn(1, :), this.dataTX(1, :), this.CSn(1))
                title('1st subcarrier')
                subplot(222), hold on, box on
                csn = find(this.CSn ~= 0, 1, 'last');
                plot_constellation(Xn(csn, :), this.dataTX(csn, :), this.CSn(csn))
                title('Last subcarrier')
                subplot(223)
                plot(abs(e(1, :)).^2)
                xlabel('Iteration')
                ylabel('MSE')
                title('Adapt. MSE for 1st subcarrier')
                subplot(224)
                plot(abs(e(end, :)).^2)
                xlabel('Iteration')
                ylabel('MSE')                
                title('Adapt. MSE for last subcarrier')
                drawnow
                
                figure(404), clf, hold on, box on
                plot(this.fc(1:length(W))/1e9, abs(W).^2, '-o')
                xlabel('Frequency (GHz)')
                ylabel('|W(f)|^2')
                title('Equalizer')
                drawnow
            end            
        end
        
        function data = demod(self, x)
            %% Demod a particular OFDM frame
            data = zeros(size(x));
            for k = 1:size(x, 1)
                if self.CSn(k) == 2
                    data(k, :) = pamdemod(x(k, :), self.CSn(k), 0, 'gray');
                elseif self.CSn(k) > 2
                    data(k, :) = qamdemod(x(k, :), self.CSn(k), 0, 'gray');
                end
            end
        end
        
        function Xn = mod(self, data)
            %% Generate a particular OFDM symbol from data
            Xn = zeros(size(data));
            for k = 1:size(data, 1)
                if self.CSn(k) == 2
                    Xn(k, :) = pammod(data(k, :), self.CSn(k), 0, 'gray');
                elseif self.CSn(k) > 2
                    Xn(k, :) = qammod(data(k, :), self.CSn(k), 0, 'gray');
                end
            end
        end
        
        function [bercount, interval] = countBER(self, Ndiscard, verbose)
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
                figure(410), box on
                stem(1:self.Nu/2, Nerr)
                xlabel('Subcarrier')
                ylabel('Number of bit errors')
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
        % 4. From that calculates this.Npre_os = Npos + Neg; (zero is not included)
        % 5. If this.Npre_os == k end simulation, otherwise increment k and repeat

%             this.Npre_os = 0;
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
            Nd = ceil(this.Nc/this.Ms);
            k = 0;
            Ncp = Inf; % cyclic prefix length
            n = -512:512;
            while Ncp > k && k < this.max_iterations
                fs_temp = this.Rs*(this.Nc + k)/Nd;
                
                tn = n/fs_temp;
                hn = interp1(t, ht, tn); % retime

                % CP based on energy
                en_frac = cumsum(hn.^2)/sum(hn.^2);
                Nn = sum(en_frac >= (1 - this.frac_incl)/2 & n < 0);
                Np = sum(en_frac <= (1 + this.frac_incl)/2 & n > 0);

                Ncp = Nn + Np;
                
                k = k + 1;
            end         

            assert(k ~= this.max_iterations, 'ofdm/cyclic_prefix: CP calculation did not converge');
            
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
                axis([1.5*this.Npre_os*[-1 1]*1/this.fs -0.5 1])
                xlabel('t (s)')
                ylabel('p(t)')
                title('Impulse response of the channel')
                drawnow
            end            
        end
        
        function ber_theory = calc_ber(self, SNRn)
            %% Calculate theoretical BER from SNR at each constellation
            ber = zeros(size(SNRn));
            for k = 1:length(SNRn)
                ber(k) = berqam(self.CSn(k), SNRn(k));
            end

            ber_theory = sum(ber.*self.bn)/sum(self.bn);
        end
        
        function [ber_est, SNRn] = estimate_ber(this, xn, Gch, noiseVar, verbose)
            %% Estimate BER from signal
            % Inputs:
            % - xn : sampled received signal
            % - Gch : channel frequency response at the subcarriers frequency
            % - noiseVar : noise variance at the subcarriers frequency

            xn = xn - mean(xn);
            Pxn = mean(abs(xn).^2);
            Pxn = Pxn*this.Nu/(this.Nc + this.Npre_os); % account for energy lost in oversampling and cyclic prefix
            
            Pref = this.var(this.Pn.*abs(Gch).^2);
            
            Pnadj = this.Pn*Pxn/Pref;
            
            % Calculate unbiased estimated SNR
            SNRn = 10*log10(Pnadj.*abs(Gch).^2./noiseVar - (this.Pn.*abs(Gch).^2./noiseVar).*(sum(noiseVar)./sum(this.Pn.*abs(Gch).^2))); 
            
            ber_est = this.calc_ber(SNRn);
            
            if exist('verbose', 'var') && verbose
                figure(301), hold on, box on
                subplot(211)
                stem(this.fc/1e9, SNRn, 'fill')
                xlabel('Frequency (GHz)')
                ylabel('Estimated SNR (dB)')
                subplot(212)
                stem(this.fc/1e9, SNRn - 10*log10(this.bn) - 10*log10(2), 'fill')
                xlabel('Frequency (GHz)')
                ylabel('Estimated Eb/N0 (dB)')    
            end
        end
    end
    
    %% Get and Set Methods
    methods 
        % OFDM sampling rate
        function fs = get.fs(this)
            if isempty(this.Npre_os)
                fs = this.Rs*(this.Nc)/this.Nu;
                warning('ofdm/fs: Cyclic prefix length was not calculated yet. Assuming no cyclic prefix.') 
            else
                fs = this.Rs*(this.Nc + this.Npre_os)/this.Nu;
            end
        end        
        
        % OFDM oversampling ratio
        function Ms = get.Ms(this)
            Ms = this.Nc/this.Nu;
        end
        
        % frequency at which subcarriers are located
        function fc = get.fc(this)
            fc = this.fs/this.Nc*(1:this.Nu/2);      
        end
        
        % Symbol rate
        function Rs = get.Rs(this)
            Rs = 2*this.Rb/log2(this.CS); 
        end
        
        % Total number of bits
        function B = get.B(this)
            B = this.Nu/2*log2(this.CS);    
        end
        
        function Npre_os = get.Npre_os(this)
            %% Cyclic prefix length
            if any(isempty([this.Npos this.Nneg]))
                Npre_os = 0;
            else
                Npre_os = this.Npos + this.Nneg;
            end
        end
        
        function Pqam = get.Pqam(this)
        %% Calculate the average power of a CS-QAM constellation with dmin = 2 (used to normalize symbol power after qammod)
            Pqam = zeros(1, this.Nu/2);
            for k = 1:this.Nu/2
                if this.CSn(k) == 0
                    continue
                else
                    Pqam(k) = mean(abs(qammod(0:this.CSn(k)-1, this.CSn(k), 0, 'gray')).^2);
                end
            end
        end  
        
        function bn = get.bn(self)
            bn = zeros(size(self.CSn));
            bn(self.CSn ~= 0) = log2(self.CSn(self.CSn ~= 0));
        end
        
        function set.powerAllocationType(self, type)
            %% Power allocation type
            if strcmpi(type, 'preemphasis') || strcmpi(type, 'palloc') || strcmpi(type, 'none')
                self.powerAllocationType = lower(type);
            else
                error('ofdm/powerAllocationType: Invalid power allocation type. It must be either preemphasis, palloc, or none.')
            end
        end
    end
end