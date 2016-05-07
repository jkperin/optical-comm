classdef AnalogLogic < AnalogOperation % inherits properties and methods from class AnalogOperation
    %% Performs XOR operation between two values. If inputs are greater and equal to 0, they are treated as '1', else they're treated as '0'. 
    %% Ideal ouptu would be -Vout or Vout. Both inputs and the output are filtered by filt
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
        
        %% XOR
        function y = xor(self, x1, x2)
            %% XOR function: output is in between -Vout and Vout. Inputs and outputs are filtered by filt.
            if self.ideal
                y = self.ideal_xor(x1, x2);
                return
            end            
            h = self.filt.h; % FIR approximation of filter
            
            % Filter inputs and multiply
            self.stateIn1 = self.shiftState(self.stateIn1, x1);
            self.stateIn2 = self.shiftState(self.stateIn2, x2);
            
            x1 = 2*(sum(self.stateIn1.*h) >= 0) - 1; % {-1, 1}
            x2 = 2*(sum(self.stateIn2.*h) >= 0) - 1; % {-1, 1}
            y = x1*x2; % this is equivalent to a XOR for values in {-1, 1}
            
            % Add thermal noise
            y = y + sqrt(self.noiseVar())*randn(1);
            
            % Filter output
            self.stateOut = self.shiftState(self.stateOut, y);
            y = self.Vout*sum(self.stateOut.*h);
        end
        
        function y = ideal_xor(self, x1, x2)
            %% Ideal XOR operation: no noise and no filtering
            x1 = 2*(x1 >= 0) - 1; % {-1, 1}
            x2 = 2*(x2 >= 0) - 1; % {-1, 1}
            y = x1*x2*self.Vout; % this is equivalent to a XOR for values in {-1, 1}
        end   
                
        function y = or(self, x1, x2)
            %% OR function: output is in between -Vout and Vout. Inputs and outputs are filtered by filt.
            if self.ideal
                y = self.ideal_or(x1, x2);
                return
            end            
            h = self.filt.h; % FIR approximation of filter
            
            % Filter inputs and multiply
            self.stateIn1 = self.shiftState(self.stateIn1, x1);
            self.stateIn2 = self.shiftState(self.stateIn2, x2);
            
            x1 = sum(self.stateIn1.*h) >= 0; % {0, 1}
            x2 = sum(self.stateIn1.*h) >= 0; % {0, 1}
            y = or(x1, x2);
            y = 2*y-1; % {-1, 1}
                        
            % Add thermal noise
            y = y + sqrt(self.noiseVar())*randn(1);
            
            % Filter output
            self.stateOut = self.shiftState(self.stateOut, y);
            y = self.Vout*sum(self.stateOut.*h);
        end
        
        function y = ideal_or(self, x1, x2)
            %% Ideal XOR operation: no noise and no filtering
            x1 = x1 >= 0; % {0, 1}
            x2 = x2 >= 0; % {0, 1}
            y = or(x1, x2);
            y = self.Vout*(2*y-1); % {-Vout, Vout}
        end 
    end
end