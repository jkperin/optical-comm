classdef DPSK < QAM
    %% Properties inherited from class QAM
%         M % DPSK order
%         Rb % bit rate
%         pulse_shape % struct containing pulse shape parameters.
%   % Dependent properties inherited from class QAM
%         Rs % symbol rate   

    properties (Dependent, GetAccess=private)
        phi0 % initial phase
    end
        
    %% methods inherited from class QAM
    % [xt, xd] = signal(self, dataTX)
    % delay = pulse_shape_grpdelay(self)
    % h = norm_filter_coefficients(~, h)
    % validate_pulse_shape(self)
    
    %% Main methods
    methods
        function obj = DPSK(M, Rb, pulse_shape)
            %% Class constructor
            if ~exist('pulse_shape', 'var')
                pulse_shape = select_pulse_shape('rect', 1);
            end
            
            obj = obj@QAM(M, Rb, pulse_shape);
        end
        
        function DPSKtable = summary(self)
            %% Generate table summarizing class values
            disp('DPSK class parameters summary:')
            rows = {'DPSK order'; 'Bit Rate'; 'Symbol rate'; 'Pulse shape'; 'Samples per symbol'};
            Variables = {'M'; 'Rb'; 'Rs'; 'pulse_shape.type'; 'pulse_shape.sps'};
            Values = {self.M; self.Rb/1e9; self.Rs/1e9; self.pulse_shape.type; self.pulse_shape.sps};
            Units = {''; 'Gbit/s'; 'Gbaud'; ''; 'S/symb'};

            DPSKtable = table(Variables, Values, Units, 'RowNames', rows);
        end
    
        function x = mod(self, dataTX)
            %% Generate symbols
            % Treats each row as an independent data stream           
            dphi = 2*pi/self.M;
            x = zeros(size(dataTX));
            for k = 1:size(dataTX, 1)
                phaseOffset = dphi*gray2bin(dataTX(k, :), 'dpsk', self.M);
                x(k, :) = sqrt(2)*exp(1j*(self.phi0 + cumsum(phaseOffset)));
            end
        end
        
        function dataRX = demod(self, x, x0)
            %% Decode symbols
            % Treats each row as an independent data stream
            % Inputs:
            % - x: symbols
            % - x0 (optional, default = phi0): previous symbol
            if not(exist('x0', 'var'))
                x0 = sqrt(2)*exp(1j*self.phi0);
            end
            
            if length(x0) == 1
                x0 = x0*ones(1, size(x, 1));
            end
              
            dphi = 2*pi/self.M;
            const = exp(1j*(dphi*gray2bin(0:self.M-1, 'dpsk', self.M)));
            distfun = @(A,B) abs(A-B).^2;
            dataRX = zeros(size(x));
            for k = 1:size(x, 1)
                phaseOffset = diff(angle([x0(k) x(k, :)]));
                dist = bsxfun(distfun, exp(1j*phaseOffset(:)), const);
                [~, databin] = min(dist, [], 2);
                dataRX(k, :) = databin.'-1;
            end
        end
        
        
        function ber = ber_freq_offset(self, SNRdB, foff)
            %% DPSK ber including frequency offset between two symbols
            % Results are based on Pawula, R. F., Rice, S. O., & Roberts, J. H. (2001). Distribution of the phase angle between two vectors perturbed by Gaussian noise II. IEEE Transactions on Vehicular Technology, 50(2), 576–583. http://doi.org/10.1109/25.923069            
            % Eqs. 9, 12, 43
            % Note: if foff = 0, then BER becomes equal to the one calculated using berawgn()
            SNR = 10.^(SNRdB/10);
            DeltaPhi = 2*pi*foff/self.Rs;
            t = linspace(-pi/2, pi/2, 1e3);
            F = @(psi, SNR) SNR*sin(DeltaPhi-psi)/(4*pi)*...
                    trapz(t, exp(-(SNR-SNR*cos(DeltaPhi-psi)*cos(t)))...
                    ./(SNR - SNR*cos(DeltaPhi-psi)*cos(t)));
            ber = zeros(size(SNR));
            for k = 1:length(SNR)
                ber(k) = 2/log2(self.M)*(F(pi, SNR(k)) - F(pi/self.M, SNR(k)));
            end
        end
    end
    
    %% Get and set methods
    methods
        %% Get and set methods
        function phi0 = get.phi0(self)
            if self.M == 4
                phi0 = pi/4; % particular case just to be compatible with QPSK constellation
            else
                phi0 = 0;
            end
        end
        
        function x = type(~)
            x = 'DPSK';
        end
       
        %% Auxiliary methods
        function validate_mod_demod(self)
            dataTX = randi([0 self.M-1], [1 2^12]);
            x = self.mod(dataTX);
            xref = sqrt(2)*exp(1j*self.phi0)*dpskmod(dataTX, self.M, 0, 'gray');

            EbN0dB = 3:15;
            EbN0 = 10.^(EbN0dB/10);
            Eb = mean(abs(self.mod(0:self.M-1)).^2)/log2(self.M);
            sig2 = Eb./EbN0;
            cropAt = 100;
            for k = 1:length(EbN0)
                yn = x + sqrt(sig2(k)/2)*(randn(size(x)) + 1j*randn(size(x)));
                yrefn = xref + sqrt(sig2(k)/2)*(randn(size(x)) + 1j*randn(size(x)));
                dataRX = self.demod(yn);
                dataRXcropped = self.demod(yn(cropAt:end), x(cropAt-1));
                dataRXref = dpskdemod(yrefn/(sqrt(2)*exp(1j*self.phi0)), self.M, 0, 'gray').';
                [~, ber(k)] = biterr(dataTX, dataRX);
                [~, berref(k)] = biterr(dataTX, dataRXref);
                [~, bercrop(k)] = biterr(dataTX(cropAt:end), dataRXcropped);
            end

            figure, hold on, box on
            plot(EbN0dB, log10(ber), '-o')
            plot(EbN0dB, log10(bercrop), '-v')
            plot(EbN0dB, log10(berref), '-s')
            plot(EbN0dB, log10(berawgn(EbN0dB, 'dpsk', self.M)), '--k')
            xlabel('E_b/N_0 (dB)')
            ylabel('log_{10}(BER)')
            legend('DSPK.mod/DSPK.mod', 'DSPK.mod/DSPK.mod with cropped sequence', 'Matlab', 'theory')
            title(['BER vs E_b/N_0 curves for ' num2str(self.M, 0) '-DPSK'])
            axis([EbN0dB(1) EbN0dB(end) -8 0])
        end    
    end   
end