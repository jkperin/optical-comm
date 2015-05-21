classdef soa
    properties
        Gain % Gain in linear units
        Fn   % noise figure (dB)
        lamb % wavelength of operation (m)
        maxGain % maximum gain (dB)
        noise_stats % shot noise statistics: {'off', 'Gaussian', 'Chi2'}, default = Gaussian
    end
    properties (Dependent)
        N0
    end
    properties (Constant)
        h = 6.62606957e-34; % Planck
        q = 1.60217657e-19; % electron charge
        c = 299792458;      % speed of light
    end
    
    methods
        %% constructor
        function obj = soa(Gain, Fn, lamb, noise_stats, maxGain)
            obj.Gain = Gain;
            obj.Fn = Fn;
            obj.lamb = lamb;
            
            if nargin == 4
                obj.noise_stats = noise_stats;
            else
                obj.noise_stats = 'Gaussian';
            end
            
            if nargin == 5
                obj.maxGain = maxGain;
            else
                obj.maxGain = Inf;
            end
        end
        
        function N0 = get.N0(obj)
            N0 = (obj.Gain - 1)*10^(obj.Fn/10)/2*(obj.h*obj.c/obj.lamb); % one-sided PSD
        end
        
        %% Set noise statistics
        function obj = set.noise_stats(obj, str)
            if strcmp(str, 'Gaussian') || strcmp(str, 'DoublyStoch') || strcmp(str, 'off')
                obj.noise_stats = str;
            else
                error('Invalid option: noise_stats should be either off, Gaussian or DoublyStoch');
            end
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
            
        
    
    