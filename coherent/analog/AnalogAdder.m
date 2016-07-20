classdef AnalogAdder < AnalogOperation % inherits properties and methods from class AnalogOperation
    methods
        function obj = AnalogAdder(filt, N0, fs)
            %% Constructor
            obj@AnalogOperation(filt, N0, fs);
        end
        
        function yf = add(self, x1, x2)
            %% Adder function: add signals and add noise. Inputs and output is filtered by filt.
            if self.ideal
                yf = self.ideal_add(x1, x2);
                return
            end            
            % Filter inputs
            [x1f, x2f] = self.filter_inputs(x1, x2);
            
            % Perform operation: comparing
            y = x1f + x2f;
            
            % Add noise
            yn = self.add_noise(y);
            
            % Filter output
            yf = self.filter_output(yn);
        end
        
        function y = ideal_add(~, x1, x2)
            %% Ideal adder operation: no noise and no filtering
            y = x1 + x2;
        end        
        
        function validate(self)
            %% Validate operation mix for the non-ideal case, where filtering is performed
            self.reset();
            N = 100;
            w = 2*pi*self.filt.fcnorm*self.fs/4;
            [~, t] = freq_time(N, self.fs);
            
            x1 = sin(w*t + pi*(2*rand(1)-1));
            x2 = sin(2*w*t + pi*(2*rand(1)-1));
            
            x1fref = filter(self.filt.num, self.filt.den, x1);
            x2fref = filter(self.filt.num, self.filt.den, x2);
            yref = self.ideal_add(x1fref, x2fref);
            ynref = self.add_noise(yref);
            yfref = filter(self.filt.num, self.filt.den, ynref);
            
            y = zeros(1, N);
            for k = 1:N
                y(k) = self.add(x1(k), x2(k));
            end
            
            figure, clf, hold on, box on
            plot(t, y)
            plot(t, yfref, '--')
            plot(t, x1+x2, ':');
            legend('this class', 'reference', 'ideal')
            self.reset();
        end        
    end
end