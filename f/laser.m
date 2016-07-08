classdef laser
    properties
        lambda % wavelength (m)
        PdBm   % Average optical power (dBm)
        RIN    % Relative intensity noise (dB/Hz)
        linewidth % Linewidth (Hz)
        freqOffset = 0;  % frequency offset with respect to wavelength lambda (Hz)
        alpha  = 0; % laser chirp
        H = @(f) ones(size(f)); % frequency response function handle
    end
    
    properties(Dependent)
        PW % power in Watts
    end      
    
    methods
        function obj = laser(lambda, PdBm, RIN, linewidth, freqOffset)
            %% Class constructor
            obj.lambda = lambda;
            obj.PdBm = PdBm;
            
            if exist('RIN', 'var')
                obj.RIN = RIN;
            end
            
            if exist('linewidth', 'var')
                obj.linewidth = linewidth;
            end
            
            if exist('freqshift', 'var')
                obj.freqOffset = freqOffset;
            end         
        end
        
        function LaserTable = summary(self)
            %% Generate table summarizing class values
            disp('-- Laser class parameters summary:')
            rows = {'Wavelength'; sprintf('Frequency offset from %.2f nm', self.lambda*1e9);...
                'Power'; 'Relative intensity noise'; 'Linewidth'};
            Variables = {'lambda'; 'freqOffset'; 'PdBm'; 'RIN'; 'linewidth'};
            Values = [self.lambda*1e9; self.freqOffset/1e9; self.PdBm; self.RIN; self.linewidth/1e3];
            Units = {'nm'; 'GHz'; 'dBm'; 'dB/Hz'; 'kHz'};

            LaserTable = table(Variables, Values, Units, 'RowNames', rows);
        end
        
        %% Get Methods
        function PW = get.PW(self)
            PW = dBm2Watt(self.PdBm);
        end
        
        %% Main Methos
        function sigma2 = varRIN(self, P, Df)
            %% RIN variance
            % Inputs
            % - P: Laser power (W)
            % - Df: One-sided bandwidth (Hz)
            sigma2 = [];
            if ~isempty(self.RIN)
                sigma2 = 10^(self.RIN/10)*P.^2*Df; 
            end
        end
        
        function sigma2 = varPN(self, fs)
            %% Variance of phase noise
            % Inputs:
            % - fs: sampling frequency
            sigma2 = [];
            if ~isempty(self.linewidth)
                sigma2 = 2*pi*self.linewidth/fs; 
            end
        end
        
        function Pout = addIntensityNoise(self, Pout, fs)
            %% Add intensity noise to signal Pin sampled at fs
            % Inputs: 
            % Pin: optical power
            % fs: sampling rate
            Pout = Pout + sqrt(self.varRIN(Pout, fs/2)).*randn(size(Pout));
        end
        
        function [Eout, phase_noise] = addPhaseNosie(self, Eout, fs)
            %% Add phase noise
            % Inputs:
            % - Ein: Electric field
            % - fs: sampling rate
            initial_phase    = pi*(2*rand(1)-1); % [-pi, pi]
            dtheta      =       [0 sqrt(self.varPN(fs))*randn(1, length(Eout)-1)]; % i.i.d.Gaussian random variables with zero mean and variance sigma_p^2
            phase_noise = initial_phase + cumsum(dtheta, 2);
            Eout = Eout.*exp(1j*phase_noise); % adds phase noise
        end
        
        function Eout = addTransientChirp(self, Ein)
            %% Add transient chirp
            Eout = Ein.*exp(1j*self.alpha/2*log(abs(Ein).^2));
        end
        
        function Eout = cw(self, sim)
            %% Generates continous wave waveform
            % Inputs:
            % sim specifies sampling rate (fs), total number of points (N),
            % and whether intensity noise (RIN) and phase noise (phase_noise)
            % are included in simulations.
            % If freqshift ~= 0, then sim must also contain time vector (t)
            
            Pout = self.PW*ones(1, sim.N);
            
            if isfield(sim, 'RIN') && sim.RIN
                Pout = self.addIntensityNoise(Pout, sim.fs);
            end
            
            Eout = sqrt(Pout);
            
            if isfield(sim, 'phase_noise') && sim.phase_noise
                Eout = self.addPhaseNosie(Eout, sim.fs);
            end
            
            if self.freqOffset ~= 0
                Eout = freqshift(Eout, sim.t, self.freqOffset);
            end            
        end    
        
        function Pout = modulate(self, x, sim)
            %% Modulate laser output based on driving signal x
            % This generates an electric signal with average power (PdBm)
            % The laser output is filtered by the frequency response H
            % If after filtering the output signal becomes negative,
            % clipping is done
            % sim specifies sampling rate (fs), total number of points (N),
            % and whether to include intensity noise (RIN) 
            % Inputs:
            %  - x: driving signal
            %  - f
            Pout = self.PW*ones(size(size(x)));
            
            Pout = real(ifft(fft(Pout).*ifftshift(self.H(sim.f))));
                                  
            if isfield(sim, 'RIN') && sim.RIN
                Pout = self.addIntensityNoise(Pout, sim.fs);
            end
            
            if self.alpha ~= 0
                Pout = self.addTransientChirp(Pout);
            end
            
            % Clip
            Pout(Pout < 0) = 0;
        end          
    end
end
