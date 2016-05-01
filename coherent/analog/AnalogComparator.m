classdef AnalogComparator < AnalogOperation % inherits properties and methods from class AnalogOperation
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
        
        function y = compare(self, x1, x2)
            %% Compare function: compare x1 and x2. If x1 >= x2, then outputs 1, else it outputs -1. Inputs and outputs are filtered by filt.
            if self.ideal
                y = self.ideal_compare(x1, x2);
                return
            end            
            h = self.filt.h; % FIR approximation of filter
            
            % Filter inputs and multiply
            self.stateIn1 = self.shiftState(self.stateIn1, x1);
            self.stateIn2 = self.shiftState(self.stateIn2, x2);
            
            y = sum(self.stateIn1.*h) >= sum(self.stateIn2.*h);
            y = (2*y - 1)*self.Vout;
            
            % Add thermal noise
            y = y + sqrt(self.noiseVar())*randn(1);
            
            % Filter output
            self.stateOut = self.shiftState(self.stateOut, y);
            y = sum(self.stateOut.*h);
        end
        
        function y = ideal_compare(self, x1, x2)
            %% Ideal adder operation: no noise and no filtering
            y = x1 >= x2;
            y = (2*y-1)*self.Vout;
        end        
    end
end