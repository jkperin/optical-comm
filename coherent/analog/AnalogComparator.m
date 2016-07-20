classdef AnalogComparator < AnalogOperation % inherits properties and methods from class AnalogOperation
    % Input1 -> Filter \
    %                   -> Compare -> Add noise -> Filter -> Output
    % Input2 -> Filter /
    % If filt is empty, component is assumed to be ideal
    properties
        Vout % output voltage
    end
    methods
        function obj = AnalogComparator(filt, N0, fs, Vout)
            %% Constructor
            obj@AnalogOperation(filt, N0, fs);
            
            if exist('Vout', 'var')
                obj.Vout = Vout;
            else
                obj.Vout = 1;
            end
        end
        
        function yf = compare(self, x1, x2)
            %% Compare function: compare x1 and x2. If x1 >= x2, then outputs Vout, else it outputs -Vout. Inputs and outputs are filtered by filt.
            if self.ideal
                yf = self.ideal_compare(x1, x2);
                return
            end   
            
            % Filter inputs
            [x1f, x2f] = self.filter_inputs(x1, x2);
            
            % Perform operation: comparing
            y = (x1f >= x2f);
            y = (2*y - 1)*self.Vout; % -Vout or Vout
            
            % Add noise
            yn = self.add_noise(y);
            
            % Filter output
            yf = self.filter_output(yn);
        end
        
        function y = ideal_compare(self, x1, x2)
            %% Ideal compare operation: no noise and no filtering
            y = x1 >= x2;
            y = (2*y-1)*self.Vout;
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
            yref = self.ideal_compare(x1fref, x2fref);
            ynref = self.add_noise(yref);
            yfref = filter(self.filt.num, self.filt.den, ynref);
            
            y = zeros(1, N);
            for k = 1:N
                y(k) = self.compare(x1(k), x2(k));
            end
            
            figure, clf, hold on, box on
            plot(t, y)
            plot(t, yfref, '--')
            plot(t, self.ideal_compare(x1, x2), ':');
            legend('this class', 'reference', 'ideal')
            self.reset();
        end  
    end
end