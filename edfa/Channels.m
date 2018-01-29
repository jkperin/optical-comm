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
        lnm
    end
    
    properties (Constant, Hidden)
        h = 6.62606957e-34; % Planck
        q = 1.60217657e-19; % electron charge
        c = 299792458;      % speed of light
    end
    
    methods
        function obj = Channels(varargin)
            %% Constructor:
            % Channels(wavelength, P, direction)
            if length(varargin) == 3 % Passed P, lamb, and direction
                obj.wavelength = varargin{1};
                obj.P = varargin{2};
                if length(obj.P) == 1
                    obj.P = obj.P*ones(size(obj.wavelength));
                end
                obj.direction = varargin{3};
            end
        end
        
        % Get/Set methods
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
        
        function lnm = get.lnm(self)
            %% wavelength in nm
            lnm = self.wavelength*1e9;
        end
        
        function E = Ephoton(self)
            %% Photon energy
            E = self.h*self.c./self.wavelength;
        end
        
        function S = sample(self, idx)
            %% Create another class Channels with a subset of channels from the orignal class
            S = Channels(self.wavelength(idx), self.P(idx), self.direction);
        end
        
        function plot(self, h)
            %% Plot channels power in dBm
            if exist('h', 'var') % handle was passed
                figure(h)
            else
                figure
            end
            
            plot(self.wavelength*1e9, self.PdBm);
            
            if not(exist('h', 'var')) 
                xlabel('Wavelength (nm)')
                ylabel('Power (dBm)')
            end
        end
    end
end