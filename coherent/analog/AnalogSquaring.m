classdef AnalogSquaring < AnalogOperation % inherits properties and methods from class AnalogOperation
    %% Output = (Input)^2
    % Input1 -> Filter -> Squanring -> Add noise -> Filter -> Output
    % If filt is empty, component is assumed to be ideal
    methods
        function obj = AnalogSquaring(filt, N0, fs)
            %% Constructor
            obj@AnalogOperation(filt, N0, fs); % calls constructor of parent class Analog Operation
        end
        
        function yf = square(self, x1)
            %% Square function: Square signal and add noise. Inputs and output is filtered by filt.
            if self.ideal
                yf = self.ideal_square(x1);
                return
            end
    
            % Filter inputs
            x1f = self.filter_inputs(x1);
            
            % Perform operation: squaring
            y = x1f.^2;
            
            % Add noise
            yn = self.add_noise(y);
            
            % Filter output
            yf = self.filter_output(yn);
        end
        
        function y = ideal_square(~, x1)
            %% Ideal squaring operation: no noise and no filtering
            y = x1.^2;
        end
        
        function validate(self)
            %% Validate operation for the non-ideal case, where filtering is performed
            self.reset();
            N = 100;
            w = 2*pi*self.filt.fcnorm*self.fs/4;
            [~, t] = freq_time(N, self.fs);
            
            x1 = sin(w*t + pi*(2*rand(1)-1));
            
            x1fref = filter(self.filt.num, self.filt.den, x1);            
            yref = self.ideal_square(x1fref);
            ynref = self.add_noise(yref);
            yfref = filter(self.filt.num, self.filt.den, ynref);
            
            y = zeros(1, N);
            for k = 1:N
                y(k) = self.square(x1(k));
            end
            
            figure, clf, hold on, box on
            plot(t, y)
            plot(t, yfref, '--')
            plot(t, x1.^2, ':');
            legend('self.mix', 'reference', 'ideal')
            self.reset();
        end
    end
end