classdef AnalogABS < AnalogOperation % inherits properties and methods from class AnalogOperation
    %% Performs ABS operation full-wave rectifier
    % Input1 -> Filter -> ABS -> Add noise -> Filter -> Output
    % If filt is empty, component is assumed to be ideal
    methods
        function obj = AnalogABS(filt, N0, fs)
            %% Constructor
            obj@AnalogOperation(filt, N0, fs); % calls parent class constructor
        end
        
        function varargout = copy(self)
            %% Deep copy of ABS. Filters states aren't copied
            for k = 1:nargout
                varargout{k} = AnalogABS(self.filt, self.N0, self.fs);
            end
        end
        
        function yf = abs(self, x1)
            %% ABS function: output is in between -Vout and Vout. Inputs and outputs are filtered by filt.
            if self.ideal
                yf = self.ideal_abs(x1);
                return
            end            
            % Filter inputs
            x1f = self.filter_inputs(x1);
            
            % Perform operation: ABS
            y = abs(x1f);
            
            % Add noise
            yn = self.add_noise(y);
            
            % Filter output
            yf = self.filter_output(yn);
        end
        
        function y = ideal_abs(~, x1)
            %% Ideal ABS operation: no noise and no filtering
            y = abs(x1);
        end   
        
        function validate(self)
            %% Validate operation for the non-ideal case, where filtering is performed
            self.reset();
            N = 100;
            w = 2*pi*self.filt.fcnorm*self.fs/4;
            [~, t] = freq_time(N, self.fs);
            
            x1 = sin(w*t + pi*(2*rand(1)-1));
            
            x1fref = filter(self.filt.num, self.filt.den, x1);            
            yref = self.ideal_abs(x1fref);
            ynref = self.add_noise(yref);
            yfref = filter(self.filt.num, self.filt.den, ynref);
            
            y = zeros(1, N);
            for k = 1:N
                y(k) = self.abs(x1(k));
            end
            
            figure, clf, hold on, box on
            plot(t, y)
            plot(t, yfref, '--')
            plot(t, abs(x1), ':');
            legend('this class', 'reference', 'ideal')
            self.reset();
        end
    end
end