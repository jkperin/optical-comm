classdef AnalogPath < AnalogOperation % inherits properties and methods from class AnalogOperation
    %% Just filters signals. This class can be used to add the same group delay as other analog operations
    % Input1 -> Filter -> Nothing -> Add noise -> Filter -> Output
    % If filt is empty, component is assumed to be ideal
    methods
        function obj = AnalogPath(filt, N0, fs)
            %% Constructor
            obj@AnalogOperation(filt, N0, fs); % calls parent class constructor
        end
        
        function varargout = copy(self)
            %% Deep copy of Path. Filters states aren't copied
            for k = 1:nargout
                varargout{k} = AnalogPath(self.filt, self.N0, self.fs);
            end
        end
        
        function yf = path(self, x1)
            %% ABS function: output is in between -Vout and Vout. Inputs and outputs are filtered by filt.
            if self.ideal
                yf = self.ideal_path(x1);
                return
            end            
            % Filter inputs
            x1f = self.filter_inputs(x1);
            
            % Perform operation: just output filtered input
            y = x1f;
            
            % Add noise
            yn = self.add_noise(y);
            
            % Filter output
            yf = self.filter_output(yn);
        end
        
        function y = ideal_path(~, x1)
            %% Ideal ABS operation: no noise and no filtering
            y = x1;
        end   
    end
end