classdef AnalogCubing < AnalogOperation % inherits properties and methods from class AnalogOperation
    %% Output = (Input)^3
    % Input1 -> Filter -> Cubing -> Add noise -> Filter -> Output
    % If filt is empty, component is assumed to be ideal
    methods
        function obj = AnalogCubing(filt, N0, fs)
            %% Constructor
            obj@AnalogOperation(filt, N0, fs); % calls constructor of parent class Analog Operation
        end
        
        function varargout = copy(self)
            %% Deep copy of Cubing. Filters states aren't copied
            for k = 1:nargout
                varargout{k} = AnalogCubing(self.filt, self.N0, self.fs);
            end
        end
        
        function yf = cube(self, x1)
            %% Cube function: cube signal and add noise. Inputs and output is filtered by filt.
            if self.ideal
                yf = self.ideal_cube(x1);
                return
            end
    
            % Filter inputs
            x1f = self.filter_inputs(x1);
            
            % Perform operation: cubing
            y = x1f.^3;
            
            % Add noise
            yn = self.add_noise(y);
            
            % Filter output
            yf = self.filter_output(yn);
        end
        
        function y = ideal_cube(~, x1)
            %% Ideal cubing operation: no noise and no filtering
            y = x1.^3;
        end
        
        function validate(self)
            %% Validate operation for the non-ideal case, where filtering is performed
            self.reset();
            N = 100;
            w = 2*pi*self.filt.fcnorm*self.fs/4;
            [~, t] = freq_time(N, self.fs);
            
            x1 = sin(w*t + pi*(2*rand(1)-1));
            
            x1fref = filter(self.filt.num, self.filt.den, x1);            
            yref = self.ideal_cube(x1fref);
            ynref = self.add_noise(yref);
            yfref = filter(self.filt.num, self.filt.den, ynref);
            
            y = zeros(1, N);
            for k = 1:N
                y(k) = self.cube(x1(k));
            end
            
            figure, clf, hold on, box on
            plot(t, y)
            plot(t, yfref, '--')
            plot(t, x1.^3, ':');
            legend('this class', 'reference', 'ideal')
            self.reset();
        end        
        
    end
end