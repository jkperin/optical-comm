classdef AnalogOperation < handle
    %% Class Analog operation: simulates an anlog operation of two inputs and generates one output
    % Input1 -> Filter \
    %                   -> Operation -> Add noise -> Filter -> Output
    % Input2 -> Filter /
    % If filt is empty, component is assumed to be ideal
    properties
        filt % filter structure created by design_filter.m function
        fs % sampling frequency of simulation
        N0 % one-sided thermal noise PSD
    end
    
    properties(Dependent)
      ideal % whether it is assumed to be ideal
    end
    
    properties
        % memory of filter in input 1, 2, and output
        input1Filter % filter of input 1
        input2Filter % filter of input 2
        outputFilter % output filter
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
            
            if not(isempty(filt))
                obj.input1Filter = ClassFilter(filt.num, filt.den, obj.fs);
                obj.input2Filter = ClassFilter(filt.num, filt.den, obj.fs);
                obj.outputFilter = ClassFilter(filt.num, filt.den, obj.fs);
            end 
        end
               
        function [x1f, x2f] = filter_inputs(self, x1, x2)               
            % Input 1
            x1f = self.input1Filter.filter(x1);
            
            % Input 2, if provided
            if exist('x2', 'var')
                x2f = self.input2Filter.filter(x2);
            else
                x2f = [];
            end
        end
        
        function yf = filter_output(self, y)
            yf = self.outputFilter.filter(y);
        end
               
        function yn = add_noise(self, y)
            %% Add white noise whose one-sided PSD is N0
            yn = y + sqrt(self.noiseVar())*randn(size(y));
        end
        
        function groupDelay = groupDelay(self)
            %% Calculates group delay in s
            if self.ideal
                groupDelay = 0;
            else
                groupDelay = 2*grpdelay(self.filt.num, self.filt.den, 1)/self.fs;
                % Note: factor of 2 appears because inputs and output are
                % filtered
            end
        end           
        
        function self = reset(self)
            %% Reset all states
            self.input1Filter.reset();
            self.input2Filter.reset();
            self.outputFilter.reset();
        end
    end
    
    methods (Access=protected)
        function sigma2 = noiseVar(self)
            sigma2 = self.N0*self.fs/2;            
        end
    end
    
    %% Get and set methods
    methods  
        function ideal = get.ideal(self)
            ideal = isempty(self.filt);
        end
    end
    
    methods
        function validate_filtering(self)
            self.reset();
            N = 100;
            x = randn(1, N);
            yref = filter(self.filt.num, self.filt.den, x);
            y = zeros(1, N);
            for t = 1:N
                y(t) = self.input1Filter.filter(x(t));
            end
            
            figure, clf, box on, hold on
            plot(x, ':')
            plot(y)
            plot(yref, '--')
            legend('Input signal', 'Output signal', 'Output signal using matlab function filter')
            self.reset();
        end          
    end
end