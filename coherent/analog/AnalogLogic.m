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
        
        function y = xnor(self, x1, x2)
            %% XOR function: output is in between -Vout and Vout. Inputs and outputs are filtered by filt.
            if self.ideal
                y = self.ideal_xnor(x1, x2);
                return
            end            
            h = self.filt.h; % FIR approximation of filter
            
            % Filter inputs and multiply
            self.stateIn1 = self.shiftState(self.stateIn1, x1);
            self.stateIn2 = self.shiftState(self.stateIn2, x2);
            
            x1 = (sum(self.stateIn1.*h) >= 0); 
            x2 = (sum(self.stateIn2.*h) >= 0); 
            y = 2*not(xor(x1, x2))-1; 
            y = self.Vout*y;
            
            % Add thermal noise
            y = y + sqrt(self.noiseVar())*randn(1);
            
            % Filter output
            self.stateOut = self.shiftState(self.stateOut, y);
            y = sum(self.stateOut.*h);
        end
        
        function y = ideal_xnor(self, x1, x2)
            %% Ideal XOR operation: no noise and no filtering
            x1 = (x1 >= 0);
            x2 = (x2 >= 0);
            y = 2*not(xor(x1, x2))-1; 
            y = self.Vout*y;
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
            y =  self.Vout*(2*y-1); % {-1, 1}
                        
            % Add thermal noise
            y = y + sqrt(self.noiseVar())*randn(1);
            
            % Filter output
            self.stateOut = self.shiftState(self.stateOut, y);
            y = sum(self.stateOut.*h);
        end
        
        function y = ideal_or(self, x1, x2)
            %% Ideal XOR operation: no noise and no filtering
            x1 = x1 >= 0; % {0, 1}
            x2 = x2 >= 0; % {0, 1}
            y = or(x1, x2);
            y = self.Vout*(2*y-1); % {-Vout, Vout}
        end         
        
        function y = and(self, x1, x2)
            %% AND function: output is in between -Vout and Vout. Inputs and outputs are filtered by filt.
            if self.ideal
                y = self.ideal_and(x1, x2);
                return
            end            
            h = self.filt.h; % FIR approximation of filter
            
            % Filter inputs and multiply
            self.stateIn1 = self.shiftState(self.stateIn1, x1);
            self.stateIn2 = self.shiftState(self.stateIn2, x2);
            
            x1 = sum(self.stateIn1.*h) >= 0; % {0, 1}
            x2 = sum(self.stateIn1.*h) >= 0; % {0, 1}
            y = and(x1, x2);
            y =  self.Vout*(2*y-1); % {-1, 1}
                        
            % Add thermal noise
            y = y + sqrt(self.noiseVar())*randn(1);
            
            % Filter output
            self.stateOut = self.shiftState(self.stateOut, y);
            y = sum(self.stateOut.*h);
        end
        
        function y = ideal_and(self, x1, x2)
            %% Ideal XOR operation: no noise and no filtering
            x1 = x1 >= 0; % {0, 1}
            x2 = x2 >= 0; % {0, 1}
            y = and(x1, x2);
            y = self.Vout*(2*y-1); % {-Vout, Vout}
        end 
        
    end
end