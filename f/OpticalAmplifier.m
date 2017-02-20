classdef OpticalAmplifier < handle
    %% Semiconductor Optical Amplifier
    properties
        Gain % Gain in linear units (same gain is assumed for both polarizations)
        Fn   % noise figure (dB)
        Wavelength % Wavelength of operation (m)
        Operation = 'ConstantOutputPower' % Amplifier operation mode = {'ConstantOutputPower', 'ConstantGain'}
        outputPower % output power in dBm. Must be set if Operation = 'ConstantOutputPower'
        maxGaindB = 30 % maximum gain (dB)
    end
    properties (Dependent)
        GaindB % Gain in dB
        Ssp % ASE one-sided PSD per real dimension
        nsp % spontaneous emission factor
    end
    properties (Constant, Hidden)
        BWref = 12.5e9; % reference bandwidth for measuring OSNR
        h = 6.62606957e-34; % Planck
        q = 1.60217657e-19; % electron charge
        c = 299792458;      % speed of light
    end
    
    methods
        function obj = OpticalAmplifier(Operation, param, Fn, Wavelength)
            %% Class constructor
            % Inputs:
            % - Opertation: either 'ConstantOutputPower' or 'ConstantGain'
            % - param: GaindB if Operation = 'ConstantGain', or outputPower
            % if Operation = 'ConstantOutputPower'
            % - Fn:  noise figure in dB
            % - Wavelength: operationl wavelength in m
            obj.set_amplifier_operation(Operation, param);
            obj.Fn = Fn;
            obj.Wavelength = Wavelength;
        end
        
        function Amptable = summary(self)
            %% Generate table summarizing class parameters
            disp('Optical amplifier class parameters summary:')
            rows = {'Operation'; 'Gain'; 'Output Power'; 'Noise Figure'; 'Operational Wavelength'};
            Variables = {'Operation'; 'GaindB'; 'outputPower'; 'Fn'; 'Wavelength'};
            Values = {self.Operation; self.GaindB; self.outputPower; self.Fn; self.Wavelength*1e9};
            Units = {''; 'dB'; 'dBm'; 'dB'; 'nm'};

            Amptable = table(Variables, Values, Units, 'RowNames', rows);
        end
        
        %% Get methods
        function Ssp = get.Ssp(obj)
            %% ASE one-sided power spectrum density per polarization
            % assumming Gain >> 1
            Ssp = (obj.Gain - 1)*10^(obj.Fn/10)/2*(obj.h*obj.c/obj.Wavelength); % Agrawal 6.1.15 3rd edition
        end
        
        function GaindB = get.GaindB(obj)
            %% Gain in dB
            GaindB = 10*log10(obj.Gain); 
        end
        
        function nsp = get.nsp(obj)
            %% Spontaneous emission factor
            nsp = 1/2*10^(obj.Fn/10);
        end 
               
        %% Set methods
        function set.GaindB(this, GdB)
            %% Gain in dB
            this.Gain = 10^(GdB/10); % set Gain, since GaindB is dependent
        end
                
        function self = set_amplifier_operation(self, operation, param)
            %% Set amplifier operation and specify Gain or output Power
            if strcmpi(operation, 'ConstantOutputPower')
                self.Operation = 'ConstantOutputPower';
                self.outputPower = param; % output power in dBm
                self.Gain = [];
            elseif strcmpi(operation, 'ConstantGain')
                self.Operation = 'ConstantGain';
                self.GaindB = param; % Gain in dB
            else
                error('OpticalAmplifier: invalid amplifier operation. Operation must be either ConstantOutputPower or ConstantGain')
            end
        end
    end
    
    %% Main Methods
    methods        
        function [varNoise, varSigSpont, varSpontSpont] = varNoiseDD(self, Pmean, Deltafele, Deltafopt, Npol)
            %% Noise variance due to ASE after direct detection. This includes signal-ASE beat noise and ASE-ASE beat noise
            % Note: this function assumes that the amplifier Gain is set,  
            % even when amplifier is operating in the constant output power mode
            % Inputs:
            % - Pmean: average signal power before amplification
            % - Deltafele: electrical noise bandwidth
            % - Deltafopt: optical noise bandwidth
            % - Npol(optical, default Npol = 2): number of noise polarizations
            
            if not(exist('Npol', 'var'))
                Npol = 2;
            end
            
            % Signal-ASE beat noise
            varSigSpont = 4*self.Gain*Pmean*self.Ssp*Deltafele;
            % ASE-ASE beat noise
            varSpontSpont = 2*Npol*self.Ssp^2*Deltafopt*Deltafele;
            
            varNoise = varSigSpont + varSpontSpont;
            % Note: Responsivity is assumed to be 1. For different
            % responsivity make sig2 = R^2*sig2
        end
        
        function adjust_gain(self, PindBm)
            %% Adjust gain to keep desire output power if amplifier is operation in constant output power mode
            if strcmpi(self.Operation, 'ConstantOutputPower')
                if isempty(self.outputPower)
                    error('OpticalAmplifier: output power must be set when amplifier Operation is ConstantOutputPower')
                end
                PoutdBm = self.outputPower;
                self.GaindB = PoutdBm - PindBm;
                if self.GaindB < 7
                    warning('OpticalAmplifier: Amplifier gain is %.2f dB. This may lead to innacuraceis since (G-1)/G is no longer approximately 1.', self.GaindB)
                elseif self.GaindB > self.maxGaindB
                    warning('OpticalAmplifier: Amplifier gain is %.2f dB. This is higher than property maxGaindB', self.GaindB)
                end
            end
        end
        
        function [Eout, OSNRdB] = amp(obj, Ein, fs)
            %% Amplification
            % - Ein = received electric field (must be a N x 1 or 2 matrix
            % depending on number of polarizations)
            % - fs = sampling frequency     
            % Ouptuts:
            % - output: amplified signal
            % - OSNRdB: OSNR in 0.1nm in dB
            
            if size(Ein, 1) == 1 % 1 pol
                Ein = [Ein; zeros(size(Ein))];
            end
            
            N = length(Ein);
            
            obj.adjust_gain(power_meter(Ein))
            
            % Ssp is the psd per polarization
            W = sqrt(obj.Ssp*fs/2)*(randn(2, N) + 1j*randn(2, N));
            % Note 1: Ssp is one-sided ASE PSD. Thus we multiply by sim.fs/2
            % Note 2: Ssp is PSD per real dimension
            
            Eout = sqrt(obj.Gain)*Ein + W;
            
            [~, Ps] = power_meter(Ein);
            obj.outputPower = power_meter(Eout);
            
            % Unbiased OSNR estimate
            OSNRdB = 10*log10(obj.Gain*Ps/(2*obj.Ssp*obj.BWref));
            
            fprintf('Optical amplifier: output power = %.2f dBm | OSNR = %.2f dB\n',...
                obj.outputPower, OSNRdB);
        end      
    end
end