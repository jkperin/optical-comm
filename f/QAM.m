classdef QAM
    properties
        M % QAM order
        Rb % bit rate
        pulse_shape % struct containing pulse shape parameters.
    end
    
    properties (Dependent)
        Rs % symbol rate
    end
       
    methods
        function obj = QAM(M, Rb, pulse_shape)
            %% Class constructor
            assert(isInteger(log2(M)), 'QAM/QAM: M must be power of 2')
            obj.M = M;
            obj.Rb = Rb;
            if exist('pulse_shape', 'var')
                obj.pulse_shape = pulse_shape;
            else
                obj.pulse_shape = select_pulse_shape('rect', 1);
            end
        end
        
        function QAMtable = summary(self)
            %% Generate table summarizing class values
            disp('QAM class parameters summary:')
            rows = {'QAM order'; 'Bit Rate'; 'Symbol rate'; 'Pulse shape'; 'Samples per symbol'};
            Variables = {'M'; 'Rb'; 'Rs'; 'pulse_shape.type'; 'pulse_shape.sps'};
            Values = {self.M; self.Rb/1e9; self.Rs/1e9; self.pulse_shape.type; self.pulse_shape.sps};
            Units = {''; 'Gbit/s'; 'Gbaud'; ''; 'S/symb'};

            QAMtable = table(Variables, Values, Units, 'RowNames', rows);
        end
    
        function x = mod(self, dataTX)
            %% Generate symbols
            % Treats each row as an independent data stream
            x = zeros(size(dataTX));
            for k = 1:size(dataTX, 1)
                x(k, :) = qammod(dataTX(k, :), self.M, 0, 'gray');
            end
        end
        
        function dataRX = demod(self, x)
            %% Decode symbols
            % Treats each row as an independent data stream
            dataRX = zeros(size(x));
            for k = 1:size(x, 1)
                dataRX(k, :) = qamdemod(x(k, :), self.M, 0, 'gray');
            end
            
        end
        
        function [xt, xd] = signal(self, dataTX)
            %% Generate QAM signal with pulse shape specified in pulse_shape.type and with pulse_shape.sps samples per symbol
            % This function treats each row as independent data streams
            % Note: group delay due to pulse shaping filtering is not removed
            % Input:
            % - dataTX = transmitted symbols from 0 to M-1
            % Outputs:
            % - xt = pulse shaped QAM signal with oversampling ratio = mpam.pulse_shape.sps (1 x length(dataTX)*mpam.pulse_shape.sps)
            % - xd = symbols at symbol rate (1 x length(dataTX))
            
            % Normalize filter taps to preserve levels amplitude           
            self.pulse_shape.h = self.norm_filter_coefficients(self.pulse_shape.h);
            % Generate data, upsample, and filter
            xd = self.mod(dataTX); % 1 x length(dataTX)
            ximp = upsample(xd.', self.pulse_shape.sps).';
            xt = filter(self.pulse_shape.h, 1, ximp, [], 2);             
        end
        
        function delay = pulse_shape_grpdelay(self)
            %% Calculate group delay of pulse shaping operation;
            delay = (length(self.pulse_shape.h)-1)/2;
        end
    end
    
    %% Get and set methods
    methods
        function Rs = get.Rs(self)
            %% Symbol rate
            Rs = self.Rb/log2(self.M);
        end
        
        function self = set.Rs(self, Rs)
            %% Set symbol rate and update bit rate
            self.Rs = Rs;
            self.Rb = Rs*log2(self.M);
        end
        
        function x = type(~)
            x = 'QAM';
        end
    end
    
    %% Auxiliary methods
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
        
        function validate_pulse_shape(self)
            Mct = self.pulse_shape.sps;
            dataTX = randi([0 self.M-1], [1 1024]);
            xt = self.signal(dataTX);
            
            figure
            subplot(211)
            eyediagram(real(xt), 2*Mct)
            title('eyediagram: I')
            subplot(212)
            eyediagram(imag(xt), 2*Mct)
            title('eyediagram: Q')
            
            figure
            self.pulse_shape.h = self.pulse_shape.h/abs(sum(self.pulse_shape.h)); % normalize to preserve levels amplitude
            [H, w] = freqz(self.pulse_shape.h, 1);
            plot(self.Rs*self.pulse_shape.sps*w/(2*pi), abs(H).^2)
            xlabel('Normalized frequency (Hz)')
            ylabel('Amplitude')
            title('Frequency response of pulse shaping filter')
        end
    end
end