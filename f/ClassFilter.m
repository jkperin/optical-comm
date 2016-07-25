classdef ClassFilter < handle
    %% Filter
    % This class can be used as a filter within a loop. It can take one 
    % value at each iteration. This is useful for simulating a PLL. When
    % used this way, the filter is implemented as a direct form 1. 
    % If a vector is passed as a parameter, then the filtering operation
    % will be performed using the filter function from Matlab (transpose
    % direct form 1), or if the flag removeGroupDelay is set, the filter
    % will be implemented in the frequency domain
    properties 
        description = {}; % cell containing order and type, if known
        num % numerator = num(1) + num(2)z^-1 + num(3)z^-2+...
        den % denominator = den(1) + den(2)z^-2 + den(3)z^-2+...
        fs % sampling frequency
        removeGroupDelay = true; % whether group delay should be removed. This is only used when a vector is passed to function filter
        % Note: group delay is calculated as the group delay at DC
    end
    
    properties(Dependent)
        filterType % FIR or IIR
        BW
    end
    
    properties
        memForward % register for forward path
        memFeedback % register for feedback path
    end
    
    methods
        function obj = ClassFilter(varargin)
            %% Constructor     
            if ischar(varargin{1})
                %% 1. filter specifications are passed: type, order, normalized cutoff frequency, fs (optional)
                type = varargin{1};
                order = varargin{2};
                fcnorm = varargin{3};
                if length(varargin) == 4
                    obj.fs = varargin{4};
                else
                    obj.fs = 1;
                end
                
                filt = design_filter(type, order, fcnorm);
                obj.num = filt.num;
                obj.den = filt.den;
                obj.description = {order, type};
            elseif isstruct(varargin{1})
                %% 2. filter struct designed with design_filter.m
                filt = varargin{1};
                type = filt.type;
                order = filt.order;
                if length(varargin) == 2
                    obj.fs = varargin{2};
                else
                    obj.fs = 1;
                end
                
                obj.num = filt.num;
                obj.den = filt.den;
                obj.description = {order, type};                
            else
                %% 3. filter is already known: passing num, den, fs (optional)
                obj.num = varargin{1};
                obj.den = varargin{2};

                if length(varargin) == 3
                    obj.fs = varargin{3};
                else
                    obj.fs = 1;
                end
            end
            
            obj.memForward = zeros(1, length(obj.num));
            obj.memFeedback = zeros(1, length(obj.den)-1);
        end
        
        function y = filter(self, x)
            %% Filter function
            if length(x) > 1 % x is the entire signal so can filter all at once
                if self.removeGroupDelay
                    delay = grpdelay(self.num, self.den, 1);
                    if strcmpi(self.filterType, 'FIR') && isInteger(delay)
                        y = filter(self.num, self.den, x);
                        if size(y, 1) > size(y, 2)
                            y = circshift(y, [delay 0]);
                        else
                           y = circshift(y, [0 delay]);
                        end
                    else
                        f = freq_time(length(x), self.fs);
                        y = ifft(fft(x).*ifftshift(self.H(f)));
                        if isreal(x) % if input is real, remove spurious imaginary component
                            y = real(y);
                        end
                    end
                else % filters using transpose direct form 1
                    y = filter(self.num, self.den, x);
                end
                return
            end

            self.memForward = self.shiftState(self.memForward, x);
            if strcmpi(self.filterType, 'FIR')
                y = sum(self.memForward.*self.num);
            else
                y = (sum(self.memForward.*self.num)- sum(self.memFeedback.*self.den(2:end)))/self.den(1);
                self.memFeedback = self.shiftState(self.memFeedback, y);
            end
        end
        
        function groupDelay = groupDelay(self)
            %% Override group delay method from parent class because only one filtering operation is performed
            groupDelay = grpdelay(self.num, self.den, 1)/self.fs;
        end
        
        function fvtool(self)
            %% Opens Matlab's fvtool for this filter
            fvtool(self.num, self.den)
        end
        
        function Hf = H(self, f)
            %% Frequency response
            Hf = freqz(self.num, self.den, f, self.fs);
            if self.removeGroupDelay
                Hf = Hf.*exp(1j*2*pi*f*self.groupDelay);
            end 
        end
        
        function self = reset(self)
            %% Reset all states
            self.memForward = zeros(size(self.memForward));
            self.memFeedback = zeros(size(self.memFeedback));
        end
        
        function newFilter = cascade(this, otherFilter)
            %% Cascade this filter with another filter
            newNum = conv(this.num, otherFilter.num);
            newDen = conv(this.den, otherFilter.den);
            newFilter = ClassFilter(newNum, newDen, this.fs);
        end
        
        function state = shiftState(~, state, x)
            %% Shift states for a new input
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
    
    %% Get and set methods
    methods
        function filterType = get.filterType(self)
            if length(self.den) == 1
                filterType = 'FIR';
            else
                filterType = 'IIR';
            end
        end
    end
end