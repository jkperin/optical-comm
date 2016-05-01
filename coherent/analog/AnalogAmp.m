classdef AnalogAmp < handle
    properties
        GaindB % Amplifier gain in dB
        filt % filter structure created by design_filter.m function
        fs % sampling frequency of simulation
        N0 % one-sided noise PSD of mixer circuitry 
    end
    properties(Dependent)
        Gain
    end
    
    properties(GetAccess=protected, Hidden)
        % memory of filter output
        stateOut
    end
    methods
        function obj = AnalogAmp(GaindB, N0, filt, fs)
            %% Constructor
            obj.GaindB = GaindB;
            
            if exist('N0', 'var')
                obj.N0 = N0;
            else 
                obj.N0 = 0;
            end
             
            obj.filt = filt;
            obj.fs = fs;
            
            obj.stateOut = zeros(size(obj.filt.h));  
        end
        
        % Get and set methods
        function Gain = get.Gain(self)
            %% Get Gain in linear units
            Gain = 10^(self.GaindB/10);
        end
        
        function set.Gain(self, G)
            %% Set Gain in linear units
            self.Gain = G;
            self.GaindB = 10*log10(G);
        end          
        
        function y = amp(self, x)
            %% Amplify 
            h = self.filt.h; % FIR approximation of filter
                      
            % Add thermal noise
            y = self.Gain*x + sqrt(self.noiseVar())*randn(1);
            
            % Filter output
            self.stateOut = self.shiftState(self.stateOut, y);
            y = sum(self.stateOut.*h);
        end
        
        function state = shiftState(~, state, x)            
            %% Shift states for a new input
            state(2:end) = state(1:end-1);
            state(1) = x;
        end            
        
        function self = reset(self)
            %% Reset all states
            self.stateOut = zeros(size(self.stateOut));
        end
    end
    
    methods (Access=protected)
        function sigma2 = noiseVar(self)
            sigma2 = self.N0*self.fs/2;            
        end
    end
end