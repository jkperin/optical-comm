classdef OSA
    %% Optical Spectrum Analyzer
    properties
        resolution % resolution in nm
    end
    
    properties(Dependent, Hidden)
        resm % resolution in m
    end
    
    properties (Constant, Hidden)
        c = 299792458;      % speed of light
        tolVec = [1e-3 1e-2 1e-1 1]; % tolerances for measuring noise floor
    end
    
    methods
        function obj = OSA(resolution)
            %% Constructor
            obj.resolution = resolution;            
        end
        
        %% Get methods       
        function resm = get.resm(self)
            resm = self.resolution*1e-9;
        end
        
        function [YdBm, ly, fy] = spectrum(self, E, lambda, f, verbose)
            %% Calculate optical spectrum
            % Inputs:
            % - x: input signal
            % - lambda: wavelength of the input signal
            % - f: frequency vector
            % - verbose (optional, default = false): whether to plot
            % spectrum at the end
            % Outputs:
            % - YdBm: power spectrum in dBm at each ly (wavelength in m) 
            % and fy (frequency in Hz) 
            
            % Calculate spectrum
            if size(E, 1) == 1 % 1 pol
                P = abs(fftshift(fft(E))).^2;
            else % 2-pol
                P = abs(fftshift(fft(E(1, :)))).^2 + abs(fftshift(fft(E(2, :)))).^2;
            end

            resolutionHz = self.resHz(lambda);
            fs = 2*max(abs(f)); % sampling rate
            P = P/(fs*length(P)/resolutionHz); 

            df = abs(f(2) - f(1));
            windowLength = round(resolutionHz/df);
            if mod(windowLength, 2) == 0 % if even
                windowLength = windowLength + 1;
            end

            Pfilt = filter(ones(1, windowLength)/windowLength, 1, P);
            Pfilt = circshift(Pfilt, [0 -(windowLength-1)/2]); % remove delay due to filter

            fc = self.c/lambda;
            fy = fc + [fliplr(0:-resolutionHz:min(f)) resolutionHz:resolutionHz:max(f)];
            Y = interp1(f, Pfilt, fy-fc);
            ly = self.c./fy;
            YdBm = 10*log10(Y/1e-3);

            if exist('verbose', 'var') && verbose
                figure(132), box on
                plot(ly*1e9, YdBm)
                xlabel('Wavelength (nm)')
                ylabel('Power (dBm)')
                title(sprintf('Optical spectrum with resolution %.2f nm', self.resolution))
            end
        end
        
        function [OSNRdB, OSNRdB01nm] = estimate_osnr(self, E, lambda, f, verbose)
            %% Estimate OSNR from optical spectrum
            % OSNR is measured assuming resolution defined in
            % self.resolution
            % IMPORTANT: measured OSNR will only match theoretical OSNR
            % when the signal bandwidth is smaller than OSA.resHz. Since,
            % this is typically not the case, the OSA should be set to a
            % course resolution (e.g., 0.3nm) so that all the signal power
            % is measured. Then, the OSNR can be converted to 0.1nm
            % resolution using the function OSA.convert_OSNR(OSNRdB, lambda)
            % IMPORTANT: this function only works correctly if the sampling 
            % rate is much larger than the symbol rate
            % Inputs:
            % - E: electric field
            % - lambda: wavelength in m
            % - f: frequency vector
            % - verbose (optinal, default = false): whether to plot the results

            % Calculate optical spectrum
            [YdBm, ly] = self.spectrum(E, lambda, f, verbose);

            [SandN, idx] = max(YdBm); % signal + noise power

            % Measure noise floor assuming different tolerances
            for k = 1:length(self.tolVec)
                tol = self.tolVec(k);
                idx1 = find(abs(diff(YdBm(idx:-1:1))) < tol); % Select points where difference was below tol
                idx2 = find(abs(diff(YdBm(idx:1:end))) < tol); % Select points where difference was below tol
                if not(isempty(idx1)) && not(isempty(idx2))
                    break
                end
            end

            if isempty(idx1) || isempty(idx2)
                warning('OSA/estimate_osnr: it was not possible to find noise floor')
                OSNRdB = [];
                return;
            end

            idx1 = idx - idx1(1);
            idx2 = idx + idx2(1);

            NdB = interp1([ly(idx1), ly(idx2)], [YdBm(idx1), YdBm(idx2)], ly(idx));

            SN = 10^(SandN/10);
            N = 10^(NdB/10);

            OSNRdB = 10*log10(SN/N - 1); % unbias OSNR estimate
            OSNRdB01nm = self.convert_OSNR(OSNRdB, lambda);

            if exist('verbose', 'var') && verbose
                figure(132), hold on
                h1 = plot(1e9*ly(idx), SandN, 'ok');
                plot(1e9*[ly(idx1(1)) ly(idx2(1))], [YdBm(idx1(1)), YdBm(idx2(1))], 'x-k')
                plot(1e9*ly(idx), NdB, 'ok')
                legend(h1, sprintf('OSNR = %.2f dB (%.2f dB @ 0.1 nm)', OSNRdB, OSNRdB01nm));
            end
        end
        
        %% Auxiliary methods
        function OSNRdB = convert_OSNR(self, OSNRdBres, lambda)
            %% Convert OSNRdB from a given resolution, to OSNRdBAdj in 0.1 nm near lambda
            OSNRdB = OSNRdBres + 10*log10(self.resHz(lambda)/12.5e9);
        end
        
        function df = resHz(self, lambda)
            %% Convert OSA resolution to Hz for wavelength = lambda
            df = self.c/(lambda - self.resm/2) - self.c/(lambda + self.resm/2);
        end            
            
    end
end