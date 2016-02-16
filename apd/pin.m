%% Positive-negative photodiode
classdef pin < apd % inherits methods and properties from apd
    %% Positive-negative photodiode inherits methods and properties from class apd
    properties
        % defined in apd.m
        % R    % responsivity
        % Id   % dark current  
        % BW   % bandwidth
    end
    methods
        function self = pin(R, Id, BW)
            if ~exist('R', 'var')
                R = 1;
            end
            
            if ~exist('BW', 'var')
                BW = Inf;
            end
            
            if ~exist('Id', 'var')
                Id = 0;
            end
            % Call constructor from APD
            % apd(GaindB, ka, BW, R, Id)
            self@apd(0, 0, BW, R, Id);
        end
        
        function output = detect(this, Ein, fs, noise_stats, N0)
            %% Detection: overwrites method from apd in order to account for 2 pol
            % doubly-stochastic noise stats case is no longer treated
            % Inputs:
            % - Ein = electric field. Could be in two pol
            % - fs = sampling frequency (Hz)
            % - noise_stats = 'gaussian', 'doubly-stochastic' (not
            % implemented), or 'no noise'
            % - N0 = one-sided PSD of thermal noise from TIA
            
            if all(size(Ein) >= 2) % two pols
                if size(Ein, 1) ~= 2
                    Ein = Ein.'; % put in 2 x N format
                end
            end
            
            Pin = sum(abs(Ein).^2, 1);
            
            switch noise_stats 
                case 'gaussian'
                    % Assuming Gaussian statistics    
                    output = this.R*this.Gain*Pin + sqrt(this.varShot(Pin, fs/2)).*randn(size(Pin));
                    
                case 'no noise'
                    % Only amplifies and filters the signal (no noise).
                    % This is used in estimating the BER
                    output = this.R*this.Gain*Pin;
                    
                otherwise 
                    error('apd>detect: Invalid Option!')
            end
            
            % Frequency
            df = fs/length(Pin);
            f = (-fs/2:df:fs/2-df).';
            
            % PD frequency response
            if ~isinf(this.BW)
                output = real(ifft(fft(output).*ifftshift(this.H(f))));
                % H has unit gain at DC
            end
                        
            % Add thermal noise if N0 was provided
            if exist('N0', 'var')
                output = output + sqrt(N0*fs/2).*randn(size(Pin)); % includes thermal noise
            end             
        end
    end
end