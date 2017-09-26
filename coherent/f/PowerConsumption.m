classdef PowerConsumption
    %% Power consumption of long-haul coherent sub-systems
    % Based on [1] Pillai, B. S. G., Sedighi, B., Guan, K., Anthapadmanabhan, 
    % N. P., Shieh, W., Hinton, K. J., & Tucker, R. S. (2014). 
    % End-to-end energy modeling and analysis of long-haul coherent 
    % transmission systems. Journal of Lightwave Technology, 32(18), 3093–3111. 
    % http://doi.org/10.1109/JLT.2014.2331086
    properties
        Rb % bit rate
        M % constellation size (default: 4)
        nb % number of bits in DSP (default: 8 bits)
        process % CMOS process in m (default: 40 nm)
        V % voltage of CMOS (default: 1 V)
    end
    
    properties (Dependent)
        % Average energy per operation in J
        Egate % gate
        Eread % write/read
        Erom % ROM
        Eadd % addition
        Emult % multiplication
    end
    
    methods
        function PowerConsumption = PowerConsumption(Rb, M, process, nb, V)
            PowerConsumption.Rb = Rb;
                        
            if exist('M', 'var')
                PowerConsumption.M = M;
            else
                PowerConsumption.M = 4;
            end
            
            if exist('process', 'var')
                PowerConsumption.process = process;
            else
                PowerConsumption.process = 40e-9; % default: 40nm CMOS
            end
            
            if exist('nb', 'var')
                PowerConsumption.nb = nb;
            else
                PowerConsumption.nb = 6;
            end
            
             if exist('V', 'var')
                PowerConsumption.V = V;
            else
                PowerConsumption.V = 0.8; % default: 40nm CMOS
            end           
        end
        
        %% Get and set methods
        function Egate = get.Egate(self)
            Egate = 1e-6*0.69*self.process*self.V^2;
        end
        
        function Eread = get.Eread(self)
            Eread = 1e-6*3.43*self.process*self.V^2;
        end
        
        function Erom = get.Erom(self)
            Erom = 1e-6*1.71*self.nb*self.process*self.V^2;
        end
        
        function Eadd = get.Eadd(self)
            Eadd = 1e-6*2.57*self.nb*self.process*self.V^2;
        end
        
        function Emult = get.Emult(self)
            Emult = 1e-6*2.57*self.nb^2*self.process*self.V^2;
        end
            
        function E = Elas(self)
            %% Laser
            % 2.5 J according to [1]
            E = 2.5/self.Rb;
        end
        
        function E = Emod(self, Vcc, vpp, Rt)
            %% Modulator
            % VCC: the driver supply voltage
            % vpp: the modulator peak to peak swing voltage
            % Rt: the termination resistance for the driver as well as the MZI
            E = 8*Vcc*vpp/(Rt*self.Rb);
        end
        
        function E = Edac(self, Fd, nd, Fs)
            %% DAC
            % Fd: DAC figure of merit (1.56 ×10?12 J/conv-step [45])
            % nd: is the DAC resolution
            % Fs: sampling frequency 
            E = 4*Fd*nd.*Fs/self.Rb;
        end
        
        function E = Epd(self, R, Vbias, Prec)
            %% Photodiode
            % R: responsivity
            % Vbias: photodioded bias voltage
            % Prec: received optical power
            E = 8*R*Vbias*Prec/self.Rb;
        end
        
        function E = Etia_agc(self)
            %% TIA-AGC
            E = 1.88/(self.Rb*log2(self.M));
        end
        
        function E = Eadc(self, Fa, nadc, Fs)
            %% ADC (for time-interleaved successive-approximation ADC)
            % Fa: ADC figure of merit
            % nadc: nominal resolution (typically ENOB + 2)
            % Fs: sampling frequency
            % Note: "High-speed ADCs such as the one studied in [65] are based 
            % on time-interleaved successive-approximation technology [62], [63]. 
            % Power consumption of such ADCs scale approximately in proportion
            % to the ADC bit resolution (unlike full flash ADCs that scale 
            % according to 2^nac) and the sampling rate [62]"
            E = 4*Fa*nadc.*Fs/self.Rb;
        end
        
        function E = Ecd(self, Ncd, ros, Nfft)
            %% CD assuming that a complex multiplication is perfomed using 3 real multiplications
            % Ncd: number of taps in CD equalizer
            % ros: oversampling ratio
            % Nfft: FFT length. Should be optimized to minimize energy
            Nn = Nfft - Ncd + 1;

            % Gate, Read, ROM, Add, Mult
            Eop = [self.Egate, self.Eread, self.Erom, self.Eadd, self.Emult];
            Nop = [self.nb*Nn*(6*log2(Nn/128)-4),...
                (4*Nn + 2*Ncd - 2)*self.nb,...
                2*Nfft*log2(Nfft)-3*Nfft+8,...
                6*Nfft*log2(Nfft)-3*Nfft+8,...
                0];
            
            E = 2*ros*sum(Nop.*Eop)/(Nn*log2(self.M)); % energy per information bit
            % Note: factor of 2 is to account for two polarizations.
        end
        
        function E = Epmd(self, Npmd, Bpmd)
            %% PMD assuming that a complex multiplication is perfomed using 3 real multiplications
            % Npmd: number of taps in CD equalizer
            % Bpmd: number of symbols per clock cycle (PMD compensation
            % block size in Table I = 256)
            
            % Gate, Read, ROM, Add, Mult
            Eop = [self.Egate, self.Eread, self.Erom, self.Eadd, self.Emult];
            Nop = [0,...
                   2*self.nb*(5*Npmd-1),...
                   4*self.nb*Bpmd,...
                   14*Npmd*(Bpmd + 1) - 2,...
                   6*Npmd*(Bpmd+1)];

           E = 2*sum(Nop.*Eop)/(Bpmd*log2(self.M)); % energy per information bit
            % Note: factor of 2 is to account for two polarizations.
        end
        
        function E = Etr(self, Btr)
            %% Timing recovery
            % Btr: number of samples input per clock cycle
            
            % Gate, Read, ROM, Add, Mult
            Eop = [self.Egate, self.Eread, self.Erom, self.Eadd, self.Emult];
            Nop = [5*self.nb*Btr,...
                self.nb*(16*Btr + 6),...
                4.5*self.nb*Btr,...
                28*Btr+1,...
                7*Btr+2];
            
            E = 2*sum(Nop.*Eop)/(Btr*log2(self.M)); % energy per information bit
            % Note: factor of 2 is to account for two polarizations.
        end
        
        function E = Ecr(self, Nf, Bcr)
            %% Carrier recovery assuming that a complex multiplication is perfomed using 3 real multiplications
            % Nf: number of symbols used in frequency estimation (32)
            % Bcr: number of samples input per clock cycle
            
            % Gate, Read, ROM, Add, Mult
            Eop = [self.Egate, self.Eread, self.Erom, self.Eadd, self.Emult];
            Nop = [8.3*self.nb*Bcr,...
                self.nb*(Bcr + 2),...
                self.nb*(10*Bcr + 2*Nf + 3 + 6*Bcr/Nf),...
                30.3*Bcr + 14*Nf - 4 - 4*Bcr/Nf,...
                7*Bcr+6*Nf];
            
            ros = 2; % oversampling ratio is assumed to be 2
            E = 2*ros*sum(Nop.*Eop)/(Bcr*log2(self.M)); % energy per information bit
            % Note: factor of 2 is to account for two polarizations.
        end
    end
end