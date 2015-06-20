classdef ofdm < handle
    properties
        Nc % number of subcarriers
        Nu % number of used subcarriers 
        CS % nominal constellation size
        Rb % bit rate
    end
    
    properties (Dependent)
        Rs % symbol rate
        Ms % OFDM oversampling ratio
        fc % subcarrier frequencies
    end
    
    properties (Dependent, GetAccess=protected)
        Pn
        dataTX
        dataRX
    end
    
    properties(Dependent, Hidden)
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
        function obj = ofdm(Nc, Nu, CS, Rs)
            obj.Nc = Nc;
            obj.Nu = Nu;
            obj.CS = CS;
            obj.Rs = Rs;          
        end
        
        %% Get methods
        % OFDM oversampling ratio
        function Ms = get.Ms(this)
            Ms = this.Nc/this.Nu;
        end
        
        % Symbol rate
        function Rs = get.Rs(this)
            Rs = 2*this.Rb/log2(this.CS); 
        end
        
        function B = get.B(this)
            B = this.Nu/2*log2(this.CS);    
        end
        
        % Calculate the average power of a CS-QAM constellation with dmin = 2 (used
        % to normalize symbol power after qammod)
        function Pqam = get.Pqam(this)
            Pqam = zeros(1, this.Nu/2);
            for k = 1:this.Nu/2
                Pqam(k) = mean(abs(qammod(0:this.CS(k)-1, this.CS(k), 0, 'gray')).^2);
            end
        end
        
        %% Power allocation and bit loading
        function power_allocation(this, tx, fiber, rx, sim)
            %% Preemphasis at the transmitter
            switch sim.type
                case 'preemphasis'
                    [this.Pn, this.CS] = preemphasis(this, tx, fiber, rx, sim);
                %% Optimal power allocation and bit loading   
                case 'palloc'
                    [this.Pn, this.CS] = palloc(this, tx, fiber, rx, sim);
                otherwise 
                    error('power_allocation: invalid option') 
            end
        end
        
        %% 
        function xt = generate_signal(this, tx, sim)
            Nsymb = sim.Nsymb;
            
            %% Generate OFDM signal at chip rate (done in DSP)
            this.dataTX = zeros(Nu/2, Nsymb);
            dataTXm = zeros(Nu/2, Nsymb);
            for kk = 1:Nu/2
                this.dataTX(kk,:) = randi([0 this.CS(kk)-1], [1 Nsymb]);              % data to be encoded (symbols columnwise)
                dataTXm(kk,:) = qammod(this.dataTX(kk,:), this.CS(kk), 0, 'gray');    % encoded QAM symbols to be modulated onto subcarriers
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
            xn = Nc*ifft(Xn, this.Nc, 1); 

            % Insert cyclic prefix
            Npre_os = this.cyclic_prefix(this, tx, rx, sim);
            xncp = [xn(end-Npre_os+1:end, :); xn]; % insert cyclic prefix

            % Parallel to serial
            xncp = reshape(xncp, 1, (Nc + Npre_os)*Nsymb); % time-domain ofdm signal w/ cyclic prefix

            %% Clipping
            % Note: both positive and negative tails are clipped
            xncpc = xncp;
            xncpc(xncp < -tx.rclip*sigtx) = -tx.rclip*sigtx;
            xncpc(xncp > tx.rclip*sigtx) = tx.rclip*sigtx;

            % Check approximations
            approx.clip_prob = 2*[qfunc(tx.rclip) mean(xncp > tx.rclip*sigtx)]; % Clipping probability

            %% Quantize at the transmitter
            if sim.quantiz
                % With quantization (clip -> quantiz -> interp)  
                % Quantiz
                yqtx = linspace(-tx.rclip*sigtx, tx.rclip*sigtx, 2^sim.ENOB);     % Quantization levels
                delta1 = abs(yqtx(2)-yqtx(1));              % tx.rclip*sigtx/(2^sim.ENOB-1);        % (2^sim.ENOB-1) because of clipping done before
                [~, xncpq] = quantiz(xncpc, yqtx(1:end-1) + delta1/2, yqtx);
            else 
                % Without quantization (clip -> interp)
                % No quantization
                xncpq = xncpc;
            end

            %% Interpolation
            % Sampling rate expansion
            xt = zeros(1, sim.Mct*Nsymb*(Npre_os + this.Nc));
            xt(1:sim.Mct:end) = xncpq;

            % Filter (interpolator + ZOH)
            xt = real(ifft(ifftshift(sim.Mct*tx.filter(sim.f/sim.fs)).*fft(xt)));
            % Note: multiply by sim.Mct because DAC filter should have gain sim.Mct
            % to compensate for 1/sim.Mct due sampling rate expansion

        end
        
        function ber = detect(this, It, rx, sim)
            % Antialiasing filter of the ADC
            It = real(ifft(fft(It).*ifftshift(rx.filter.H(sim.f/sim.fs)))); % signal + noise
            % Note: we will treat noise separately as it makes it easier to estimate
            % SNR (etc). If all the operations that follow are linear this should make
            % no difference

            %% Back to discrete time
            % Resample signal once per chip        
            yncp = It(1:sim.Mct:end);       % yn with cyclic prefix without noise  
            
            %% Quantize at the receiver
            if sim.quantiz    
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
                [~, yncpq] = quantiz(yncpc, yqrx(1:end-1) + delta2/2, yqrx);
            else
                % No quantization
                yncpq = yncp - mean(yncp);
            end

            %%
            % reshape into matrix form
            yncpq = reshape(yncpq, Nc + Npre_os, Nsymb);      % signal + noise

            % Remove cyclic prefix
            yn = circshift(yncpq(Npos_os+1:end-Nneg_os, :), -Nneg_os);  % signal + noise

            % Demodulate symbols 1/N*DFT()
            Yn = fft(yn, Nc, 1)/Nc;             % signal + noise      

            % used subcarrier amplitudes (complex conjugate subcarriers are ignored)
            dataRXm = Yn(1 + (1:Nu/2), this.Ndis+1:end-this.Ndis);            

            % % AGC caculates what should be the scaling factor to normalize all
            % constellations so that dmin = 2 (value expected by qamdemod)
            K = 1 - qfunc(tx.rclip) - qfunc(rx.rclip);
            
            % Remove group delay of modulator frequency response           
            Hmod = tx.kappa*tx.modulator.H(this.fc);
            Hmod = Hmod.*exp(1j*2*pi*this.fc*tx.modulator.grpdelay);

            % The group delay is removed from the frequency response of the ADC and DAC
            Gdac = tx.filter.H(this.fc/sim.fs);                    
            Gadc = rx.filter.H(this.fc/sim.fs);   

            % Frequency response of the channel at the subcarriers
            Gch = K*Gdac.*tx.kappa.*Hmod.*Hfiber.*rx.R.*Gadc;         
            
            AGCn = sqrt(this.Pqam)./(K*sqrt(this.Pn).*Gch);
            % Note: Factor of K appears due to clipping

            this.dataTX = this.dataTX(:, this.Ndis+1:end-this.Ndis);
            this.dataRX = zeros(this.Nu/2, sim.Nsymb-2*this.Ndis);
            numerr = zeros(1, Nu/2);
            bn = log2(this.CS);
            for kk = 1:Nu/2
            %     scatterplot(dataRXm(kk,:))
                dataRXm(kk,:) = AGCn(kk)*dataRXm(kk,:);
            %     scatterplot(dataRXm(kk,:))
                this.dataRX(kk, :) = qamdemod(dataRXm(kk,:), this.CS(kk), 0, 'gray');         % encoded QAM symbols to be modulated onto subcarriers
                numerr(kk) = biterr(this.dataTX(kk, :), this.dataRX(kk,:), bn(kk));
            end

            % average BER
            ber.count(navg) = sum(numerr)/(sim.Nsymb*sum(bn));

            % 95% confidence intervals for the counted BER
            [~, interval] = berconfint(numerr, sim.Nsymb*bn);
            ber.interval(navg, :) = [min(interval(:,1)), max(interval(:,2))];

            % Estimate BER from the SNR at each subcarrier
            % Note that to estimate the BER we use SNRn which is calculated from the
            % data measured in the simulation. If we had used SNRnest (calculated from
            % the theoretical values) the results would match the target BER perfectly
%             berest = zeros(size(SNRn));
%             for kk = 1:length(SNRn)
%                 berest(kk) = berqam(this.CS(kk), SNRn(kk));
%             end
% 
%             ber.est(navg) = sum(berest.*bn)/sum(bn);
        end
               
        %% Calculate the cyclic prefix length after oversampling 
        % The channel is assumed to be a 2nd-order filter with damping ratio = 1.

        % Algorithm:
        % 1. Assumes a cyclic prefix length k
        % 2. From k, it calculates the new sampling rate fs
        % 3. Calculate the number of samples at both sides (Nneg and Npos) that 
        % contains the desired fraction of energy
        % 4. From that calculates Npre_os = Npos + Neg; (zero is not included)
        % 5. If Npre_os == k end simulation, otherwise increment k and repeat

        function [Npre_os, Nneg, Npos] = cyclic_prefix(this, tx, rx, sim)
        
            Mct = sim.Mct;
            Ntot = 1024;

            Nd = ceil(this.N/this.ros);
            
            gdac = impz(tx.filter.den, tx.filter.num);
            gadc = impz(rx.filter.den, rx.filter.num);

            k = 0;
            Npre_os = Inf;
            while Npre_os > k && k < this.max_iterations
                fs = this.Rs*(this.Nc + k)/Nd;
                fsct = Mct*fs;
                dt = 1/fsct;

                % Group delay of modulator in samples
                modulator_grpdelay = tx.modulator.grpdelay*fsct;

                % Channel impulse response
                tct = (0:Ntot-1)*dt;
                hct = tx.hl(tct);
                
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
                Nneg = sum(en_frac >= (1 - this.frac_incl)/2 & t < 0);
                Npos = sum(en_frac <= (1 + this.frac_incl)/2 & t > 0);

                % Number of samples necessary to attain desired fraction of energy
                Npre_os = Nneg + Npos;

                k = k + 1;
            end

            assert(k ~= max_it, 'CP calculation did not converge');
        end
    end
end