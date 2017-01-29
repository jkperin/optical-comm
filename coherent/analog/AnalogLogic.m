classdef AnalogLogic < AnalogOperation % inherits properties and methods from class AnalogOperation
    %% Performs logic operation between two values. If inputs are greater and equal to 0, they are treated as '1', else they're treated as '0'. 
    %% Ideal output would be -Vout or Vout. Both inputs and the output are filtered by filt
    properties
        Vout % output voltage
    end
    methods
        function obj = AnalogLogic(filt, N0, fs, Vout)
            %% Constructor
            obj@AnalogOperation(filt, N0, fs); % calls parent class constructor
            
            if exist('Vout', 'var')
                obj.Vout = Vout;
            else
                obj.Vout = 1;
            end
        end
        
        function varargout = copy(self)
            %% Deep copy of Logic. Filters states aren't copied
            for k = 1:nargout
                varargout{k} = AnalogLogic(self.filt, self.N0, self.fs);
                varargout{k}.Vout = self.Vout;
            end
        end
        
        function yf = xor(self, x1, x2)
            %% XOR function: output is in between -Vout and Vout. Inputs and outputs are filtered by filt.
            if self.ideal
                yf = self.ideal_xnor(x1, x2);
                return
            end                       
            % Filter inputs
            [x1f, x2f] = self.filter_inputs(x1, x2);
            
            % Perform operation: XOR
            x1d = (x1f >= 0); 
            x2d = (x2f >= 0); 
            y = 2*xor(x1d, x2d)-1; 
            y = self.Vout*y;
            
            % Add noise
            yn = self.add_noise(y);
            
            % Filter output
            yf = self.filter_output(yn);
        end
        
        function y = ideal_xor(self, x1, x2)
            %% Ideal XOR operation: no noise and no filtering
            x1 = (x1 >= 0);
            x2 = (x2 >= 0);
            y = 2*xor(x1, x2)-1; 
            y = self.Vout*y;
        end  
        
        function yf = xnor(self, x1, x2)
            %% XNOR function: output is in between -Vout and Vout. Inputs and outputs are filtered by filt.
            if self.ideal
                yf = self.ideal_xnor(x1, x2);
                return
            end                       
            % Filter inputs
            [x1f, x2f] = self.filter_inputs(x1, x2);
            
            % Perform operation: XNOR
            x1d = (x1f >= 0); 
            x2d = (x2f >= 0); 
            y = 2*not(xor(x1d, x2d))-1; 
            y = self.Vout*y;
            
            % Add noise
            yn = self.add_noise(y);
            
            % Filter output
            yf = self.filter_output(yn);
        end
        
        function y = ideal_xnor(self, x1, x2)
            %% Ideal XNOR operation: no noise and no filtering
            x1 = (x1 >= 0);
            x2 = (x2 >= 0);
            y = 2*not(xor(x1, x2))-1; 
            y = self.Vout*y;
        end   
                
        function yf = or(self, x1, x2)
            %% OR function: output is in between -Vout and Vout. Inputs and outputs are filtered by filt.
            if self.ideal
                yf = self.ideal_or(x1, x2);
                return
            end            
            % Filter inputs
            [x1f, x2f] = self.filter_inputs(x1, x2);
            
            % Perform operation: OR
            x1d = (x1f >= 0); 
            x2d = (x2f >= 0); 
            y = 2*or(x1d, x2d)-1; 
            y = self.Vout*y;
            
            % Add noise
            yn = self.add_noise(y);
            
            % Filter output
            yf = self.filter_output(yn);
        end
        
        function y = ideal_or(self, x1, x2)
            %% Ideal XOR operation: no noise and no filtering
            x1 = x1 >= 0; % {0, 1}
            x2 = x2 >= 0; % {0, 1}
            y = or(x1, x2);
            y = self.Vout*(2*y-1); % {-Vout, Vout}
        end         
        
        function yf = and(self, x1, x2)
            %% AND function: output is in between -Vout and Vout. Inputs and outputs are filtered by filt.
            if self.ideal
                yf = self.ideal_and(x1, x2);
                return
            end            
            % Filter inputs
            [x1f, x2f] = self.filter_inputs(x1, x2);
            
            % Perform operation: AND
            x1d = (x1f >= 0); 
            x2d = (x2f >= 0); 
            y = 2*and(x1d, x2d)-1; 
            y = self.Vout*y;
            
            % Add noise
            yn = self.add_noise(y);
            
            % Filter output
            yf = self.filter_output(yn);
        end
        
        function y = ideal_and(self, x1, x2)
            %% Ideal XOR operation: no noise and no filtering
            x1 = x1 >= 0; % {0, 1}
            x2 = x2 >= 0; % {0, 1}
            y = and(x1, x2);
            y = self.Vout*(2*y-1); % {-Vout, Vout}
        end 
        
        
        function validate_xor(self)
            %% Validate operation for the non-ideal case, where filtering is performed
            self.reset();
            N = 100;
            w = 2*pi*self.filt.fcnorm*self.fs/4;
            [~, t] = freq_time(N, self.fs);
            
            x1 = sin(w*t + pi*(2*rand(1)-1));
            x2 = sin(2*w*t + pi*(2*rand(1)-1));
            
            x1fref = filter(self.filt.num, self.filt.den, x1);
            x2fref = filter(self.filt.num, self.filt.den, x2);
            yref = self.ideal_xor(x1fref, x2fref);
            ynref = self.add_noise(yref);
            yfref = filter(self.filt.num, self.filt.den, ynref);
            
            y = zeros(1, N);
            for k = 1:N
                y(k) = self.xor(x1(k), x2(k));
            end
            
            figure, clf, hold on, box on
            plot(t, y)
            plot(t, yfref, '--')
            plot(t, self.ideal_xor(x1, x2), ':');
            legend('this class', 'reference', 'ideal')
            self.reset();
        end
    end
end