classdef soa < handle
    properties
        Gain % Gain in linear units
        Fn   % noise figure (dB)
        lamb % wavelength of operation (m)
        maxGaindB % maximum gain (dB)
    end
    properties (Dependent)
        GaindB
        N0 % psd of ASE noise (complex Gaussian) per polarization 
    end
    properties (Constant, GetAccess=protected)
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
        
        %% Get methods
        function N0 = get.N0(obj)
            % assumming Gain >> 1
            N0 = (obj.Gain - 1)*10^(obj.Fn/10)/2*(obj.h*obj.c/obj.lamb); % one-sided PSD
        end
        
        function GaindB = get.GaindB(obj)
            GaindB = 10*log10(obj.Gain); % one-sided PSD
        end
        
        %% Set methods
        function set.GaindB(this, GdB)
            this.Gain = 10^(GdB/10); % set Gain, since GaindB is dependent
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
        
        %% Optimize SOA gain
        function optimize_gain(this, mpam, tx, fiber, rx, sim)
            disp('Optimizing SOA gain...')

            if strcmp(mpam.level_spacing, 'uniform')
                % Optmize gain for uniform spacing: find Gain that minimizes the required
                % average power (Prec) to achieve a certain target SER.
                [Gsoa_opt, ~, exitflag] = fminbnd(@(Gsoa) fzero(@(PtxdBm) calc_soa_ber(PtxdBm, Gsoa, mpam, tx, fiber, this, rx, sim) - sim.BERtarget, -20), 1, min(10^(this.maxGaindB/10), 1000));    
                
            elseif strcmp(mpam.level_spacing, 'nonuniform')
                warning('optimize_gain for nonuniform spacing not implemented yet')
            end
            
            if exitflag ~= 1
                warning('SOA gain optimization did not converge (exitflag = %d)\n', exitflag)
            end 

            if ~isnan(Gsoa_opt) && ~isinf(Gsoa_opt) && Gsoa_opt >= 1
                this.Gain = Gsoa_opt;
            else
                this.Gain = originalGain;
                warning('SOA gain was not changed')
            end

            function ber = calc_soa_ber(PtxdBm, Gsoa, mpam, tx, fiber, soa, rx, sim)
                % Set power level
                tx.Ptx = 1e-3*10^(PtxdBm/10);

                % Set APD gain
                soa.Gain = Gsoa; % linear units
                
                % Uniform level spacing
                mpam.a = (0:2:2*(mpam.M-1)).';
                mpam.b = (1:2:(2*(mpam.M-1)-1)).';

                % Estimated BER using KLSE Fourier and saddlepoint approximation of
                % tail probabilities
                ber = ber_soa_klse_fourier(rx.U_fourier, rx.D_fourier, rx.Fmax_fourier, mpam, tx, fiber, soa, rx, sim);
            end
        end
     
    end
end
            
        
    
    