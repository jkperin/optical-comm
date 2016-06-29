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
    end

    properties (Dependent)
        fs % sampling rate
        Rs % symbol rate
        Ms % OFDM oversampling ratio
        fc % subcarrier frequencies
    end
    
    properties (GetAccess=private)
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
        Ndis = 16;
    end
    
    methods
        %% Constructor
        function obj = ofdm(Nc, Nu, CS, Rb, Npos, Nneg)
            obj.Nc = Nc;
            obj.Nu = Nu;
            obj.CS = CS;
            obj.Rb = Rb;     
            if nargin >= 5
                obj.Npos = Npos;
                obj.Nneg = Nneg;
            end
        end
               
        %% Power allocation and bit loading
        function power_allocation(this, tx, fiber, rx, sim)
            %% Preemphasis at the transmitter
            switch sim.type
                case 'preemphasis'
                    [this.Pn, this.CSn] = preemphasis(this, tx, fiber, rx, sim);
                %% Optimal power allocation and bit loading   
                case 'palloc'
                    [this.Pn, this.CSn] = palloc(this, tx, fiber, rx, sim);
                otherwise 
                    error('power_allocation: invalid option') 
            end
            
            if isfield(sim, 'verbose') && sim.verbose
                figure
                subplot(211), box on, grid on
                stem(1:length(this.CSn), log2(this.CSn))
                xlabel('Subcarrier')
                ylabel('log_2(CS)')
                
                subplot(212), box on, grid on
                stem(1:length(this.CSn), this.Pn, 'fill')
                xlabel('Subcarrier')
                ylabel('Power (W)');
            end
        end
        
        %% Generate OFDM signal
        function xt = generate_signal(this, tx, sim)
            Nsymb = sim.Nsymb;
            
            %% Generate OFDM signal at chip rate (done in DSP)
            this.dataTX = zeros(this.Nu/2, Nsymb);
            dataTXm = zeros(this.Nu/2, Nsymb);
            for kk = 1:this.Nu/2
                this.dataTX(kk,:) = randi([0 this.CSn(kk)-1], [1 Nsymb]);              % data to be encoded (symbols columnwise)
                dataTXm(kk,:) = qammod(this.dataTX(kk,:), this.CSn(kk), 0, 'gray');    % encoded QAM symbols to be modulated onto subcarriers
                dataTXm(kk,:) = sqrt(this.Pn(kk))*dataTXm(kk,:)/sqrt(this.Pqam(kk));  % scale constellation so that Pn is the power at nth subcarrier (E(|Xn|^2) = Pn)
            end          
            
            % zero-pad and ensure Hermitian symmetry
            % -f(N) ... -f(1) f(0) f(1) ... f(N-1) => (0 ... 0, X(-n)*, 0, X(n), 0 ... 0)
            Xn = ifftshift([zeros((this.Nc-this.Nu)/2, Nsymb); ...     % zero-pad neg. frequencies
                          flipud(conj(dataTXm));...            % data*(-n) (Nu) 
                          zeros(1, sim.Nsymb); ...                % 0 at f == 0 (1)
                          dataTXm; ...                         % data(n)  (Nu)
                          zeros((this.Nc-this.Nu)/2 - 1, Nsymb)], 1);   % zero-pad pos. frequencies

            % Perform ifft (Nc-IDFT) columnwise
            xn = this.Nc*ifft(Xn, this.Nc, 1); 

            % Insert cyclic prefix           
            xncp = [xn(end-this.Npre_os+1:end, :); xn]; % insert cyclic prefix

            % Parallel to serial
            xncp = reshape(xncp, 1, (this.Nc + this.Npre_os)*Nsymb); % time-domain ofdm signal w/ cyclic prefix

            %% Clipping
            sigtx = sqrt(2*sum(this.Pn));
            % Note: both positive and negative tails are clipped
            xncpc = xncp;
            xncpc(xncp < -tx.rclip*sigtx) = -tx.rclip*sigtx;
            xncpc(xncp > tx.rclip*sigtx) = tx.rclip*sigtx;

            %% Quantize at the transmitter
            if isfield(sim, 'quantiz') && sim.quantiz
                % With quantization (clip -> quantiz -> interp)  
                % Quantiz
                yqtx = linspace(-tx.rclip*sigtx, tx.rclip*sigtx, 2^sim.ENOB);     % Quantization levels
                delta1 = abs(yqtx(2)-yqtx(1));              % tx.rclip*sigtx/(2^sim.ENOB-1);        % (2^sim.ENOB-1) because of clipping done before
                [~, xncpq, this.varQtx] = quantiz(xncpc, yqtx(1:end-1) + delta1/2, yqtx); 
            else 
                % Without quantization (clip -> interp)
                % No quantization
                this.varQtx = 0;
                xncpq = xncpc;
            end

            %% Interpolation
            % Sampling rate expansion
            xt = zeros(1, sim.Mct*Nsymb*(this.Npre_os + this.Nc));
            xt(1:sim.Mct:end) = xncpq;

            % Filter (interpolator + ZOH)
            xt = real(ifft(ifftshift(sim.Mct*tx.filter.H(sim.f/sim.fs)).*fft(xt.')));
            % Note: multiply by sim.Mct because DAC filter should have gain sim.Mct
            % to compensate for 1/sim.Mct due sampling rate expansion

        end
        
        %% Detect OFDM signal and calculate BER
        function ber = detect(this, It, Gch, rx, sim)
            % Antialiasing filter of the ADC
            It = real(ifft(fft(It).*ifftshift(rx.filter.H(sim.f/sim.fs)))); % signal + noise
            % Note: we will treat noise separately as it makes it easier to estimate
            % SNR (etc). If all the operations that follow are linear this should make
            % no difference

            %% Back to discrete time
            % Resample signal once per chip        
            yncp = It(1:sim.Mct:end);       % yn with cyclic prefix without noise  
            
            %% Quantize at the receiver
            if isfield(sim, 'quantiz') && sim.quantiz    
                % signal std
                sigrx = std(It);
                
                % remove dc bias
                yncp = yncp - mean(yncp);

                % quantization levels
                yqrx = linspace(-rx.rclip*sigrx, rx.rclip*sigrx, 2^sim.ENOB); % Quantization levels
                delta2 = abs(yqrx(2) - yqrx(1)); 

                % clip before quantization
                yncpc = yncp;
                yncpc(yncp < yqrx(1)) = yqrx(1);
                yncpc(yncp > yqrx(end)) = yqrx(end);

                % Quantize
                [~, yncpq, this.varQrx] = quantiz(yncpc, yqrx(1:end-1) + delta2/2, yqrx);
            else
                % No quantization
                this.varQrx = 0;
                yncpq = yncp - mean(yncp);
            end

            %%
            % reshape into matrix form
            yncpq = reshape(yncpq, this.Nc + this.Npre_os, sim.Nsymb);      % signal + noise

            % Remove cyclic prefix
            yn = circshift(yncpq(this.Npos+1:end-this.Nneg, :), -this.Nneg);  % signal + noise

            % Demodulate symbols 1/N*DFT()
            Yn = fft(yn, this.Nc, 1)/this.Nc;             % signal + noise      

            % used subcarrier amplitudes (complex conjugate subcarriers are ignored)
            dataRXm = Yn(1 + (1:this.Nu/2), this.Ndis+1:end-this.Ndis);            
          
            AGCn = sqrt(this.Pqam)./(sqrt(this.Pn).*Gch);
            % Note: Factor of K appears due to clipping

            this.dataRX = zeros(this.Nu/2, sim.Nsymb-2*this.Ndis);
            numerr = zeros(1, this.Nu/2);
            bn = log2(this.CSn);
            for kk = 1:this.Nu/2
%                 scatterplot(dataRXm(kk,:))
                dataRXm(kk,:) = AGCn(kk)*dataRXm(kk,:);
%                 scatterplot(dataRXm(kk,:))
                this.dataRX(kk, :) = qamdemod(dataRXm(kk,:), this.CSn(kk), 0, 'gray');         % encoded QAM symbols to be modulated onto subcarriers
                numerr(kk) = biterr(this.dataTX(kk, this.Ndis+1:end-this.Ndis), this.dataRX(kk,:), bn(kk));
            end

            % average BER
            ber.count = sum(numerr)/(sim.Nsymb*sum(bn));

            % 95% confidence intervals for the counted BER
            [~, interval] = berconfint(numerr, sim.Nsymb*bn);
            ber.interval = [min(interval(:,1)), max(interval(:,2))];     
        end
               
        %% Calculate the cyclic prefix length after oversampling 
        % The channel is assumed to be a 2nd-order filter with damping ratio = 1.

        % Algorithm:
        % 1. Assumes a cyclic prefix length k
        % 2. From k, it calculates the new sampling rate fs
        % 3. Calculate the number of samples at both sides (Nneg and Npos) that 
        % contains the desired fraction of energy
        % 4. From that calculates this.Npre_os = Npos + Neg; (zero is not included)
        % 5. If this.Npre_os == k end simulation, otherwise increment k and repeat

        function cyclic_prefix(this, tx, rx, sim)
        
            Mct = sim.Mct;
            Ntot = 1024;

            Nd = ceil(this.Nc/this.Ms);
            
            gdac = impz(tx.filter.num, tx.filter.den);
            gadc = impz(rx.filter.num, rx.filter.den);

            k = 0;
            Ncp = Inf; % cyclic prefix length
            while Ncp > k && k < this.max_iterations
                fs_temp = this.Rs*(this.Nc + k)/Nd;
                fsct = Mct*fs_temp;
                dt = 1/fsct;

                % Group delay of modulator in samples
                modulator_grpdelay = tx.modulator.grpdelay*fsct;

                % Channel impulse response
                tct = (0:Ntot-1)*dt;
                hct = tx.modulator.h(tct);
                
                % Total impulse response in continuous time
                pct = conv(conv(gdac, hct, 'full'), gadc, 'full');
                pct = pct/max(pct);

                % Remove group delay due to filters so that impulse response is
                % centered at zero.
                tct = (0:length(pct)-1)*dt;
                tct = tct - ceil(tx.filter.grpdelay + rx.filter.grpdelay + modulator_grpdelay)*dt;

                n0 = ceil(tx.filter.grpdelay + rx.filter.grpdelay + modulator_grpdelay) + 1; % new zero after removing group delay

                % Sampling at the chip rate
                t = tct([fliplr(n0:-Mct:1), n0+Mct:Mct:end]);
                p = pct([fliplr(n0:-Mct:1), n0+Mct:Mct:end]);

                % CP based on energy
                en_frac = cumsum(p.^2)/sum(p.^2);
                this.Nneg = sum(en_frac >= (1 - this.frac_incl)/2 & t < 0);
                this.Npos = sum(en_frac <= (1 + this.frac_incl)/2 & t > 0);

                Ncp = this.Nneg + this.Npos;
                
                k = k + 1;
            end         

            assert(k ~= this.max_iterations, 'CP calculation did not converge');
            
            if isfield(sim, 'verbose') && sim.verbose
                figure, grid on, hold on, box on
                plot(tct, pct)
                stem(t, p, 'fill')
                plot(-this.Nneg/this.fs*[1 1], [-1 1], 'k')
                plot(this.Npos/this.fs*[1 1], [-1 1], 'k')
                text(-(this.Nneg)/this.fs, 0.7, sprintf('N_{neg} = %d', this.Nneg))
                text((this.Npos+0.1)/this.fs, 0.7, sprintf('N_{pos} = %d', this.Npos))
                axis([1.5*this.Npre_os*[-1 1]*1/this.fs -0.5 1])
                xlabel('t (s)')
                ylabel('p(t)')
                title('Impulse response of the channel (DAC * Laser * ADC)')
            end            
        end
        
        % Estimate BER from the SNR at each subcarrier
        % A = signal power at each subcarrier
        % varthermal = thermal noise variance of each carrier after ADC filter
        % varshot = shot noise variance of each carrier after ADC filter
        % varrin  = RIN noise variance of each carrier after ADC filter
        % Gch = total channel frequency response
        function ber_est = estimate_ber(this, A, varthermal, varshot, varrin, Gch)
            SNRn = 10*log10(this.Nc*(A^2)*this.Pn.*abs(Gch).^2./...
                (varthermal + varshot + varrin + this.varQtx*abs(Gch).^2 + this.varQrx)); % Measured SNR
            
            berest = zeros(size(SNRn));
            for kk = 1:length(SNRn)
                berest(kk) = berqam(this.CSn(kk), SNRn(kk));
            end

            bn = log2(this.CSn);
            ber_est = sum(berest.*bn)/sum(bn);    
        end
    end
    
    %% Get Methods
    methods 
        % OFDM sampling rate
        function fs = get.fs(this)
            try
                fs = this.Rs*(this.Nc + this.Npre_os)/this.Nu;
            catch e
                disp('Cyclic prefix was not calculated yet.')
                disp(e.message)
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
        
        % Cyclic prefix length
        function Npre_os = get.Npre_os(this)
            Npre_os = this.Npos + this.Nneg;
        end
        
        % Calculate the average power of a CS-QAM constellation with dmin = 2 (used
        % to normalize symbol power after qammod)
        function Pqam = get.Pqam(this)
            Pqam = zeros(1, this.Nu/2);
            for k = 1:this.Nu/2
                Pqam(k) = mean(abs(qammod(0:this.CSn(k)-1, this.CSn(k), 0, 'gray')).^2);
            end
        end   
    end
end