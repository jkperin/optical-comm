classdef ClassFilter < handle
    properties 
        ideal=true % whether it is assumed to be ideal
        num % numerator = num(1) + num(2)z^-1 + num(3)z^-2+...
        den % denominator = den(1) + den(2)z^-2 + den(3)z^-2+...
        fs % sampling frequency
    end
    
    properties(SetAccess=private)
        implementationType % FIR or IIR
        memForward % register for forward path
        memFeedback % register for feedback path
    end
    
    methods
        function obj = ClassFilter(num, den, fs)
            %% Constructor     
            obj.num = num;
            obj.den = den;
            obj.fs = fs;
            
            obj.memForward = zeros(1, length(num));
            obj.memFeedback = zeros(1, length(den)-1);
            if isempty(obj.memFeedback)
                obj.implementationType = 'FIR';
            else
                obj.implementationType = 'IIR';
            end
        end
        
        function y = filter(self, x)
            %% Filter function
            if length(x) > 1 % x is the entire signal so can filter all at once
                y = filter(self.filt.num, self.filt.den, x);
                return
            end

            self.memForward = self.shiftState(self.memForward, x);
            if strcmpi(self.implementationType, 'FIR')
                y = sum(self.memForward.*self.num);
            else
                y = (sum(self.memForward.*self.num)- sum(self.memFeedback.*self.den(2:end)))/self.den(1);
                self.memFeedback = self.shiftState(self.memFeedback, y);
            end
        end
        
        function groupDelay = groupDelay(self)
            %% Override group delay method from parent class because only one filtering operation is performed
            if self.ideal
                groupDelay = 0;
            else
                groupDelay = grpdelay(self.num, self.den, 1)/self.fs;
            end
        end
        
        function self = reset(self)
            %% Reset all states
            self.memForward = zeros(size(self.memForward));
            self.memFeedback = zeros(size(self.memFeedback));
        end
        
        
        function state = shiftState(~, state, x)
            % Shift states for a new input
            state(2:end) = state(1:end-1);
            state(1) = x;
        end       
        
        function outputSignal = validate(self, inputSignal)
            %% Check if filter implementation is working properly
            if not(exist('inputSignal', 'var'))
                inputSignal = randn(1, 1000); % generates random signal for testing
            end
            
            self.reset();
            outputSignal = zeros(size(inputSignal));
            self = self.reset();
            for t = 1:length(inputSignal)
                outputSignal(t) = self.filter(inputSignal(t));
            end
            
            figure, hold on, box on
            plot(inputSignal)
            plot(outputSignal, 'r')
            plot(filter(self.num, self.den, inputSignal), '--k');
            legend('Input signal', 'Output signal', 'Output signal using matlab function filter')
            self.reset();
        end  
    end
end