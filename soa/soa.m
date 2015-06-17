classdef soa
    properties
        Gain % Gain in linear units
        Fn   % noise figure (dB)
        lamb % wavelength of operation (m)
        maxGaindB % maximum gain (dB)
    end
    properties (Dependent)
        GaindB
        N0
    end
    properties (Constant)
        h = 6.62606957e-34; % Planck
        q = 1.60217657e-19; % electron charge
        c = 299792458;      % speed of light
    end
    
    methods
        %% constructor
        function obj = soa(GaindB, Fn, lamb, maxGaindB)
            obj.Gain = 10^(GaindB/10);
            obj.Fn = Fn;
            obj.lamb = lamb;
                        
            if nargin == 4
                obj.maxGaindB = maxGaindB;
            else
                obj.maxGaindB = Inf;
            end
        end
        
        function N0 = get.N0(obj)
            % assumming Gain >> 1
            N0 = (obj.Gain - 1)*10^(obj.Fn/10)/2*(obj.h*obj.c/obj.lamb); % one-sided PSD
        end
        
        function GaindB = get.GaindB(obj)
            GaindB = 10*log10(obj.Gain); % one-sided PSD
        end
               
        %% Amplification
        % Ein = received electric field in one pol; fs = sampling frequency
        function [output, w] = amp(obj, Ein, fs)
            % noise            
            N = length(Ein);
            
            % N0 is the psd per polarization
            w = sqrt(1/2*obj.N0*fs/2)*(randn(N, 1) + 1j*randn(N, 1));
            
            output = Ein*sqrt(obj.Gain) + w;         
        end
            
    end
end
            
        
    
    