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
            w = sqrt(1/2*obj.N0*fs)*(randn(N, 1) + 1j*randn(N, 1));
            % Note: even though soa.N0 is single-sided PSD we don't multiply by
            % sim.fs/2 because this is a band-pass process
            
            output = Ein*sqrt(obj.Gain) + w;         
        end
        
        %% Optimize SOA gain
        function optimize_gain(this, mpam, tx, fiber, rx, sim)
            disp('Optimizing SOA gain...')

            originalGain = this.Gain;
            
            % Optmize gain: find Gain that minimizes the required
            % average power (Prec) to achieve a certain target SER.
            [Gsoa_opt, ~, exitflag1] = fminbnd(@(Gsoa) calc_opt_PtxdBm(Gsoa, mpam, tx, fiber, this, rx, sim),...
                1, min(10^(this.maxGaindB/10), 1000));    
            
            if exitflag1 ~= 1
                warning('soa>optimize_gain: SOA gain optimization did not converge (exitflag1 = %d)\n', exitflag1)
            end 

            if ~isnan(Gsoa_opt) && ~isinf(Gsoa_opt) && Gsoa_opt >= 1
                this.Gain = Gsoa_opt;
            else
                this.Gain = originalGain;
                warning('soa>optimize_gain: SOA gain was not changed')
            end
          
            function PtxdBm_opt = calc_opt_PtxdBm(Gsoa, mpam, tx, fiber, soa, rx, sim)
                
                ber = zeros(size(tx.PtxdBm));
                for k = 1:length(tx.PtxdBm)
                    % Set power level
                    tx.Ptx = 1e-3*10^(tx.PtxdBm(k)/10);
                    
                    % Set SOA gain
                    soa.Gain = Gsoa; % linear units
                    
                    switch mpam.level_spacing
                        case 'uniform'
                            % Uniform level spacing
                            mpam.a = (0:2:2*(mpam.M-1)).';
                            mpam.b = (1:2:(2*(mpam.M-1)-1)).';
                        case 'nonuniform'
                            if isfield(mpam, 'level_spacing_with_gaussian_approx') && mpam.level_spacing_with_gaussian_approx
                                % Optimize level spacing using Gaussian approximation
                                Deltaf = rx.elefilt.noisebw(sim.fs)/2; % electric filter one-sided noise bandwidth
                                varTherm = rx.N0*Deltaf; % variance of thermal noise

                                Deltafopt = rx.optfilt.noisebw(sim.fs); % optical filter two-sided noise bandwidth
                                % function to calculate noise std
                                calc_noise_std = @(Plevel) sqrt(varTherm + 2*Plevel*soa.N0*Deltaf + 2*soa.N0^2*Deltafopt*Deltaf*(1-1/(2*sim.M)));
                                % Note: Plevel corresponds to the level after SOA amplification.
                                % Therefore, the soa.Gain doesn't appear in the second term because
                                % it's already included in the value of Plevel.
                                % Note: second term corresponds to sig-sp beat noise, and third term
                                % corresponds to sp-sp beat noise with noise in one polarization. Change to
                                % 2 to 4 in third term to simulate noise in two pols.
                                
                                [mpam.a, mpam.b] = level_spacing_optm_gauss_approx(mpam.M, sim.BERtarget, tx.rexdB, calc_noise_std, sim.verbose);
                            else
                                [mpam.a, mpam.b] = level_spacing_optm(mpam, tx, soa, rx, sim);
                            end
                        otherwise
                            error('soa>optimize_gain: mpam.level_spacing invalid option')
                    end

                    % Estimated BER using KLSE Fourier and saddlepoint approximation of
                    % tail probabilities
                    ber(k) = ber_soa_klse_fourier(rx.U_fourier, rx.D_fourier, rx.Fmax_fourier, mpam, tx, fiber, soa, rx, sim);
                end
                
                PtxdBm_opt = interp1(log10(ber), tx.PtxdBm, log10(sim.BERtarget), 'spline');
            end
                    
                
%                 [PtxdBm_opt, fval, exitflag] = fzero(@(PtxdBm) calc_soa_ber(PtxdBm, Gsoa, mpam, tx, fiber, soa, rx, sim) - sim.BERtarget, -20);
%             
%                 if exitflag ~= 1 % most likely is in an error floor
%                     exitflag
%                     fval
%                 end        
            
            function ber = calc_soa_ber(PtxdBm, Gsoa, mpam, tx, fiber, soa, rx, sim)
                % Set power level
                tx.Ptx = 1e-3*10^(PtxdBm/10);

                % Set APD gain
                soa.Gain = Gsoa; % linear units
                
                switch mpam.level_spacing
                    case 'uniform'
                        % Uniform level spacing
                        mpam.a = (0:2:2*(mpam.M-1)).';
                        mpam.b = (1:2:(2*(mpam.M-1)-1)).';
                    case 'nonuniform'
                        [mpam.a, mpam.b] = level_spacing_optm(mpam, tx, soa, rx, sim);
                    otherwise
                        error('soa>optimize_gain: mpam.level_spacing invalid option')
                end

                % Estimated BER using KLSE Fourier and saddlepoint approximation of
                % tail probabilities
                ber = ber_soa_klse_fourier(rx.U_fourier, rx.D_fourier, rx.Fmax_fourier, mpam, tx, fiber, soa, rx, sim);
            end
        end
     
    end
end
            
        
    
    