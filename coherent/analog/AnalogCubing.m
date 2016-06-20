classdef AnalogCubing < AnalogOperation % inherits properties and methods from class AnalogOperation
    %% Output = (Input)^3
    methods
        function obj = AnalogCubing(filt, N0, fs)
            %% Constructor
            obj@AnalogOperation(filt, N0, fs); % calls constructor of parent class Analog Operation
        end
        
        function y = cube(self, x1)
            %% Cube function: cube signal and add noise. Inputs and output is filtered by filt.
            if self.ideal
                y = self.ideal_cube(x1);
                return
            end
            
            h = self.filt.h; % FIR approximation of filter
            
            % Filter inputs and multiply
            self.stateIn1 = self.shiftState(self.stateIn1, x1);
            
            y = sum(self.stateIn1.*h)^3;
            
            % Add thermal noise
            y = y + sqrt(self.noiseVar())*randn(1);
            
            % Filter output
            self.stateOut = self.shiftState(self.stateOut, y);
            y = sum(self.stateOut.*h);
        end
        
        function y = ideal_cube(~, x1)
            %% Ideal mixer operation: no noise and no filtering
            y = x1^3;
        end
    end
end