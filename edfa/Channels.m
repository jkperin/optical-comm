classdef Channels
    properties
        wavelength % wavelength
        P % power
        direction % direction of propgation: {'forward', 'backward'}
    end
    
    properties(Dependent)
        N % number of channels
        u % direction of propagation u = 1, if direction = 'forward', u = -1, if direction = 'backward'
        PdBm % Power in dBm
    end
    
    methods
        function obj = Channels(varargin)
            %% Constructor
            if length(varargin) == 3 % Passed P, lamb, and direction
                obj.wavelength = varargin{1};
                obj.P = varargin{2};
                if length(obj.P) == 1
                    obj.P = obj.P*ones(size(obj.wavelength));
                end
                obj.direction = varargin{3};
            end
        end
        
        % Get methods
        function N = get.N(self)
            %% number of channels
            N = length(self.P);
        end
        
        function u = get.u(self)
            %% Direction of propagation
            u = ones(1, self.N);
            if strcmpi(self.direction, 'backward')
                u = -u;
            end
        end
        
        function PdBm = get.PdBm(self)
            %% Power in dBm
            PdBm = 10*log10(self.P/1e-3);
        end
    end
end