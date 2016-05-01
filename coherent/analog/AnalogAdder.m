classdef AnalogAdder < AnalogOperation % inherits properties and methods from class AnalogOperation
    methods
        function obj = AnalogAdder(filt, N0, fs)
            %% Constructor
            obj@AnalogOperation(filt, N0, fs);
        end
        
        function y = add(self, x1, x2)
            %% Adder function: add signals and add noise. Inputs and output is filtered by filt.
            if self.ideal
                y = self.ideal_add(x1, x2);
                return
            end            
            h = self.filt.h; % FIR approximation of filter
            
            % Filter inputs and multiply
            self.stateIn1 = self.shiftState(self.stateIn1, x1);
            self.stateIn2 = self.shiftState(self.stateIn2, x2);
            
            y = sum(self.stateIn1.*h) + sum(self.stateIn2.*h);
            
            % Add thermal noise
            y = y + sqrt(self.noiseVar())*randn(1);
            
            % Filter output
            self.stateOut = self.shiftState(self.stateOut, y);
            y = sum(self.stateOut.*h);
        end
        
        function y = ideal_add(~, x1, x2)
            %% Ideal adder operation: no noise and no filtering
            y = x1 + x2;
        end        
    end
end