classdef AnalogABS < AnalogOperation % inherits properties and methods from class AnalogOperation
    %% Performs ABS operation full-wave rectifier
    methods
        function obj = AnalogABS(filt, N0, fs)
            %% Constructor
            obj@AnalogOperation(filt, N0, fs); % calls parent class constructor
        end
        
        function y = abs(self, x1)
            %% ABS function: output is in between -Vout and Vout. Inputs and outputs are filtered by filt.
            if self.ideal
                y = self.ideal_abs(x1);
                return
            end            
            h = self.filt.h; % FIR approximation of filter
                      
            % Add thermal noise
            y = abs(x1) + sqrt(self.noiseVar())*randn(1);
            
            % Filter output
            self.stateOut = self.shiftState(self.stateOut, y);
            y = sum(self.stateOut.*h);
        end
        
        function y = ideal_abs(~, x1)
            %% Ideal ABS operation: no noise and no filtering
            y = abs(x1);
        end        
        
        function groupDelay = groupDelay(self)
            %% Override group delay method from parent class because only one filtering operation is performed
            if self.ideal
                groupDelay = 0;
            else
                groupDelay = grpdelay(self.filt.h, 1, 1)/self.fs;
            end
        end    
        
    end
end