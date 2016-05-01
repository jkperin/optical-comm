classdef AnalogOperation < handle
    %% Class Analog operation: simulates an anlog operation of two inputs and generates one output. This is the parent class for adder and mix
    properties
        ideal=false % whether it is assumed to be ideal
        filt % filter structure created by design_filter.m function
        fs % sampling frequency of simulation
        N0 % one-sided thermal noise PSD 
        groupDelay % group delay of operation in seconds
    end
    properties(GetAccess=protected, Hidden)
        memoryLength = 0 % memory length of state filters
        % memory of filter in input 1, 2, and output
        stateIn1 
        stateIn2
        stateOut
    end
    methods
        function obj = AnalogOperation(filt, N0, fs)
            %% Constructor
            % Assign filt, N0, and fs and initialize states
            % If filt==[], then assumes ideal simulation
            
            obj.filt = filt;
            obj.fs = fs;
            
            if isempty(N0)
               obj.N0 = 0;
            else 
               obj.N0 = N0;
            end
                        
            if isempty(filt)
                obj.ideal = true;
            end
            
            obj.stateIn1 = zeros(1, obj.memoryLength);  
            obj.stateIn2 = zeros(1, obj.memoryLength);
            obj.stateOut = zeros(1, obj.memoryLength);
        end
        
        %% Get methods
        function groupDelay = get.groupDelay(self)
            %% Calculates group delay in s
            if self.ideal
                groupDelay = 0;
            else
                groupDelay = 2*grpdelay(self.filt.h, 1, 1)/self.fs;
                % Note: factor of 2 appears because inputs and output are
                % filtered
            end

        end    
        
        function memoryLength = get.memoryLength(self)
            if isempty(self.filt)
                memoryLength = 0;
            else
                memoryLength = length(self.filt.h);
            end
        end
        
        function state = shiftState(~, state, x)
            % Shift states for a new input
            state(2:end) = state(1:end-1);
            state(1) = x;
        end            
        
        function self = reset(self)
            % Reset all states
            self.stateIn1 = zeros(1, self.memoryLength);
            self.stateIn2 = zeros(1, self.memoryLength);
            self.stateOut = zeros(1, self.memoryLength);
        end
    end
    
    methods (Access=protected)
        function sigma2 = noiseVar(self)
            sigma2 = self.N0*self.fs/2;            
        end
    end
end