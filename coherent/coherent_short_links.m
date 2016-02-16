% Ngi_CoOFDM.m     Milad Sharif
% No-Gaurd Interval Co-OFDM with M subcarriers per sampling and one
% modulator for group of M subcarriers.
% PMD, PDL, higher-order modulation, laser phase noise, fiber nonlinearity
% are not included 
% thermal noise, LO shot noise are neglected.
close all
%clear all
addpath([pwd '\plots_scripts'])
addpath([pwd '\DSP'])
%% ====== Simulation Parameter =======                                                      
DumpData = 0;                                                               % Dump data for adap. equ. optimization - make sure to be zero if you dont want to dump the data (takes a while to dumb all of te data) 
Save = 1;                                                                   % Save the final result        
PMD  = 0;
PhaseNoise   = 0;
SNR = 12;                                                                   % Optimize parameters and plot figures @SNR = SNR
FreqError = 1;
%% ============ CONSTANTS ============
lambda = 1550e-9;                                                           % wavelength (m)
c      = 3.00e8;                                                            % speed of light (m/s)
h      = 6.63e-34;                                                          % Planck's constant (J*s)
q      = 1.602e-19;                                                         % electronic charge (J)
nu     = c/lambda;                                                          % optical frequency (Hz)

%% =============== PBRS ===============
PRBS.N    = 14;                                                             % Number of bits in LFSR
%PRBS.tabs = [13,12,10,9];                                                  % LFSR taps 

%% =========== Constellation ==========
Const.type = 'QAM';
Const.M    = 4;

%% ======== Co-OFDM parameteres =======
Rs     = 12.5e9;                                                            % Symbol rate per carrier (symbol/s)
Ts     = 1/Rs;                                                              % Symbol interval per carrier (s)
M      = 1;                                                                 % Number of carriers per sampling
Ncar   = 24;                                                                % Number of carriers in superchannel (must be a multiple of M)
f1     = -325e9;                                                            % Frequency of first carrier (Hz)
fcars  = f1 + ((1:Ncar)-1)*Rs;                                              % Frequencies of all carriers (Hz)
K      = 1:Ncar/M;                                                              % Superchannels indexes   
nSC    = 4;                                    
% ----------- Frequency Error ---------
pattern = 2;                                                                % 1 : df on the desired channel
                                                                            % 2 : df on the adjacent channels
error_ratio = 20;
df      = error_ratio*Rs*1e-2; 

df_LO   = 0*Rs*1e-4;                                                        % Freq. error of LO laser 
%% ===================================
% For modulation, propagation and detection of carriers -------------------
Nsymb = 2^PRBS.N-1;                                                         % Length of PRBS
Nstep = 24;                                                                 % Time steps per symbol interval (mulitple of oversamp)
Ntot  = Nstep*Nsymb;                                                        % total number of time steps in waveforms                
t     = (0:Ntot-1)*Ts/Nstep;                                                % Continuous time vector
omega = 2*pi/(Ts/Nstep)/Ntot*[0:(Ntot/2-1), (-Ntot/2):(-1)];                % continuous freq. vector
% For digital signal processing -------------------------------------------
oversamp      = 6;
Ntotdsp  = oversamp*Nsymb;                                                  % total number of samples
tdsp     = Ts/oversamp*(0:Ntotdsp-1);                                       % time vector
omegadsp = 2*pi/(Ts/oversamp)/Ntotdsp*[0:(Ntotdsp/2-1), (-Ntotdsp/2):(-1)]; % freq. vector
deltat   = Ts/Nstep;

%% =========== Mod. ==================
ratio      = 1;%.96;                                                            
n_r        = 2.15;                                                          
n_m        = n_r*ratio; 
Mod.d_12   = (n_m -n_r)/c;                                                  % Walkoff parameter between microwave and optical fields
%Mod.a      = .01*100*sqrt(abs(omega)/2/pi*1e-9);                            % Microwave attenuation coefficient (electrode loss)
Mod.a      = .005*100*sqrt(abs(omega)/2/pi*1e-9);                            % Microwave attenuation coefficient (electrode loss)
Mod.L      = .05;                                                           % Modulator length
Mod.Hel    = (1-exp(Mod.L*(-Mod.a + 1j*Mod.d_12*omega)))./...               % Freq. response of the optical modulator limited by elec. loss and velocity mismatch
                (Mod.a-1j*Mod.d_12*omega)/Mod.L;
Mod.Hel(1) = 1;
%% =========== Transmitter ========== 
Tx.const = Const;                                                           % Constellation
Tx.mod   = Mod;                                                             % Optical modulator
Tx.Pulse = 'NRZ';                                                           % Transmitter pulse shape ('NRZ')
Tx.BW    = 1.2*Rs;                                                          % Bandwidth at -3 dB of mod. driver + mod. (Hz)
Tx.Type  = 'Bessel';                                                        % Response of mod. driver + mod. (bes or but)
Tx.Ord   = 5;                                                               % Order of transmitter composite response (1-5)
Tx.Vpi   = 4.5;                                                             % Modulator switching voltage (V)
Tx.Delx  = 0;                                                               % Delay of x pol. in pol. mux. (symbol intervals)
Tx.Dely  = 0;                                                               % Delay of y pol. in pol. mux. (symbol intervals)

%% ============ Phase Noise ==========
independent= 1;
num_lasers = 2+2*nSC;
PN.d_nu    = 200e3;                                                         % beat linewidth
PN.sigma_p = sqrt(2*pi*PN.d_nu*Ts/Nstep);                                   % sqrt of variance of phase difference over one symbol
% Generate random received phase following a Wiener process ---------------
Theta_0    = pi*(2*rand(num_lasers,1)-1);
PN.nu      = PN.sigma_p*randn(num_lasers,Ntot-1);                           % i.i.d.Gaussian random variables with zero mean and variance sigma_p^2
PN.Theta = [Theta_0 repmat(Theta_0,1,Ntot-1)+cumsum(PN.nu,2)];
if (independent==0)
    PN.Theta = [repmat(PN.Theta(1,:),1+2*nSC,1);PN.Theta(end,:)];
end
if (PhaseNoise==0)
    PN.Theta = 0*PN.Theta;
end
%% ============= Add/Drop ============
Add.ILdB  = -7.0;                                                           % Add WSS average insertion loss (dB)
Add.T     = 10^(Add.ILdB/10);
Add.H     = sqrt(Add.T)*ones(size(omega));
Drop.ILdB = -7.0;                                                           % Drop WSS average insertion loss (dB)
Drop.T    = 10^(Drop.ILdB/10);
Drop.H    = sqrt(Drop.T)*ones(size(omega));

%% =========== Mux and demux ============
MUX.FWHM = 50;                                                              % Mux and demux FWHM (GHz)
MUX.ord  = 6;                                                               % Mux and demux super-Gaussian order
MUX.H    = exp(-log(2)/2*(omega/pi/MUX.FWHM/1e9).^(2*MUX.ord));             % super-Gaussian transfer function

%% ====== Fiber (excluding PMD) ========
% Loss --------------------------------------------------------------------
Fiber.alphadBkm = -0.25;                                                    % SMF loss (dBm/km)
% Chromatic Dispersion ----------------------------------------------------
Fiber.Dpsnmkm = 17;                                                         % SMF CD (ps/(nm*km))
Dsmm          = 1e-6*Fiber.Dpsnmkm;           
Fiber.beta2   = -lambda^2/(2*pi*c)*Dsmm;                                    % SMF beta2
% Polarization-Mode Dispersion --------------------------------------------
Fiber.taumeanps = 4.5;                                                       % Mean DGD (ps)
Fiber.Nsect     = 10;                                                       % number of sections per span used to generate PMD

if (PMD==1)
    jm_filename = sprintf('jm%dps%dsec',ceil(Fiber.taumeanps),Fiber.Nsect);           % file name for Jones matrix
    if(~exist([jm_filename '.mat'],'file'))
        display('*** Creating Jones Matrix ...');
        Fiber.M = jonesMatrix(Fiber.taumeanps,Fiber.Nsect,omega);
        display('Done.')
        %save(jm_filename,'M');
    else
        display('*** Loading Jones Matrix ...');
        load(jm_filename);
        Fiber.M = PMD_M;
    end
else
    Fiber.M = zeros(2,2,length(omega));
    Fiber.M(1,1,:) = ones(size(omega));
    Fiber.M(2,2,:) = ones(size(omega));
end


% -------------------------------------------------------------------------
Fiber.Nspan   = 25; 
Fiber.Lkm     = 80;                                                         % SMF span length (km)
Fiber.Tf      = 10^(Fiber.alphadBkm*Fiber.Lkm/10);                          % Power transmittance of fiber segment
Lseg          = 1e3*Fiber.Lkm;                                              % Fiber segment length (m)
Fiber.H       = sqrt(Fiber.Tf)*exp(-1i*Fiber.beta2/2*Lseg*omega.^2);        % Freq resp of fiber segment (excluding PMD)

%% =========== Amplifier ==============
% EDFAs
NFdB      = 5.0;                                                             % EDFA noise figure (dB)
NF        = 10^(NFdB/10);               
Amp1.nsp  = NF/2;                                                            % EDFA inversion factor
Amp1.GdB  = abs(Add.ILdB+Drop.ILdB);                                         % Second EDFA gain (dB)
Amp2.nsp  = NF/2;                                                            % EDFA inversion factor
Amp2.GdB  = abs(Fiber.alphadBkm*Fiber.Lkm);                                  % Second EDFA gain (dB)

%% =========== Noise loading =============
Att.min    = -37;                                                            % Attenuation at minimum SNR (dB)
Att.max    = -20;                                                            % Attenuation at maximum SNR (dB)      

Att.stepdB = 1;                                                             % Attentuation step (dB)
Att.att = Att.min:Att.stepdB:Att.max;                                       % Attenuation for noise loading

%% =============== Receiver ===============
% Local oscillator --------------------------------------------------------
Rx.PLOdBm = 13;                                                             % Total local oscillator power (dBm)
% polarization splitting --------------------------------------------------
Rx.PolSplit.sig  = 'PBS';                                                   % pbs: polarization beamsplitter
Rx.PolSplit.LO   = 'PBS';                                                   % 3dB: 3-dB coupler     
Rx.PolSplit.Rext = 30;                                                      % PBS extinction ratio (dB), default = 30
Rx.PolSplit.R3dB = 1/2;                                                     % power splitting ratio of nominally 3-dB coupler (system performance should be insensitive to this parameter)
% 90-degree hybrid, same parameter for two polarizations ------------------
Rx.Hybrid.fS = 0.5;                                                         % power splitting ratio for signal coupler (W/W), default = 0.5
Rx.Hybrid.fL = 0.5;                                                         % power splitting ratio for LO coupler (W/W), default = 0.5
Rx.Hybrid.fI = 0.5;                                                         % power splitting ratio for in-phase coupler (W/W), default = 0.5
Rx.Hybrid.fQ = 0.5;                                                         % power splitting ratio for quadrature coupler (W/W), default = 0.5
Rx.Hybrid.tauIps = 0;                                                       % delay in in-phase branch (ps), default = 0
Rx.Hybrid.tauQps = 0;                                                       % delay in quadrature branch (ps), default = 0
Rx.Hybrid.phiI01deg = 90;                                                   % d.c. phase shift in I branch of pol. 1 (degrees), default = 90
Rx.Hybrid.phiQ01deg = 0;                                                    % d.c. phase shift in Q branch of pol. 1 (degrees), default = 0
Rx.Hybrid.phiI02deg = 90;                                                   % d.c. phase shift in I branch of pol. 2 (degrees), default = 90
Rx.Hybrid.phiQ02deg = 0;                                                      % d.c. phase shift in Q branch of pol. 2 (degrees), default = 0
% Photodiode --------------------------------------------------------------
RX.PD.Eta   = 0.8;                                                          % Photodiode quantum efficiency
Rx.PD.R     = RX.PD.Eta*q*lambda/h/c;                                       % Photodiode responsivity
Rx.LPF.BW   = 2.4*Rs;                                                       % Bandwidth at -3 dB of rec. filt.
Rx.LPF.Type = 'Butter';                                                     % Response of receiver filter
Rx.LPF.Ord  = 5;                                                            % Order of receiver filter (1-5)

%% ============== Sampler ==============
%oversamp      = 6; 

%% ============= Quantizer ==============
Quantz.Range = 8;                                                           % Range between max and min levels (V)
Quantz.Bits  = 16; 

%% ========== Digital Filters ===========
% CD compensation ---------------------------------------------------------
wintype = 'rectwin';                                                        % Window type. Options include: bartlett, chebwin, hamming, hanning, rectwin (see MATLAB "help window"). 
Ltot    = 1e3*Fiber.Lkm*Fiber.Nspan;                                                    % total SMF length (m)
eval(['win = window(@' wintype ', length(omegadsp));']);
Hcdcomp = exp(1i*Fiber.beta2/2*Ltot*omegadsp.^2).*win.';                    % Freq resp of CD compensating digital filter
% Carrier separation ------------------------------------------------------
HT2nd = 0.5*(exp(1i*omegadsp*Ts/4)+exp(-1i*omegadsp*Ts/4));                 % delay by T/2 and add (shifted to give zero group delay)
HT4nd = 0.5*(exp(1i*omegadsp*Ts/8)+exp(-1i*omegadsp*Ts/8));                 % delay by T/4 and add (shifted to give zero group delay)
Hcarsep = HT2nd.*HT2nd.*HT4nd;                                              % overall carrier separation filter
HSINC = sinc(Ts*omegadsp/2/pi);

Hcarsep = HSINC;
%% ========= Adaptive Equalizer ==========
Equalization = 'LMS';                                                       % 'CMA'    : Constant Modulus Algorithm
                                                                            % 'MMA'    : Multi-Modulus Algorithm 
                                                                            % 'RDE'    : Radius-Directed Equalization
                                                                            % 'LMS'    : Decision-Directed Least-Mean-Squared Algorithm  
                                                                            % 'CMA-DD' : CMA for pre-convergence and switch to DD after the converge
                                                                            
% Constant Modulus Algorithm ----------------------------------------------
AdEq.CMA.mu    = 0.002;
AdEq.CMA.LFilt = 6;
%getOptVal;

% Multi Modulus Algorithm -------------------------------------------------
AdEq.MMA.mu    = 0.01;
AdEq.MMA.LFilt = 3;

% Radius-Directed Equalization --------------------------------------------
AdEq.RDE.mu    = 0.01;
AdEq.RDE.LFilt = 3;

% Decision-Directed Least Mean Squared Algorithm --------------------------
AdEq.LMS.LFilt    = 10; 
AdEq.LMS.muInit   = 0.05;                                          
AdEq.LMS.muFactor = .6;
AdEq.LMS.muInt    = 1000;                                           
AdEq.LMS.muShifts = 10;
%getOptVal;

AdEq.LMS.muVec = AdEq.LMS.muInit*AdEq.LMS.muFactor.^(min((AdEq.LMS.muShifts+1),...      % vector of step size
    ceil((1:Nsymb)/AdEq.LMS.muInt))-1);     
% Training mode vs. decision-directed mode ----------------------------
AdEq.LMS.startdd = round(AdEq.LMS.muShifts/2)*AdEq.LMS.muInt;               % symbol index at which to stop using training symbols and change to decision-directed mode
AdEq.LMS.decwtvec = (1:Nsymb)<AdEq.LMS.startdd;                             % vector of decision weights (1 when use training symbols, 0 when use decision symbols)
AdEq.LMS.startmeas = AdEq.LMS.muShifts*AdEq.LMS.muInt;                               % symbol index at which to start measuring error, Q, BER

% Decision-Directed with CMA preconvergence
AdEq.CMA_DD.CMA_mu     = 0.01;
AdEq.CMA_DD.LFilt      = 7;
AdEq.CMA_DD.hardSwitch = 3000;
AdEq.CMA_DD.muInit     = 0.4; % .001~.01 %                                          
AdEq.CMA_DD.muFactor   = .6; % .2~.7 %
AdEq.CMA_DD.muInt      = 500; % 500~1000 %                                          
AdEq.CMA_DD.muShifts   = 14; % 6~20 %
%getOptVal;


AdEq.CMA_DD.muInit
AdEq.CMA_DD.muVec   = [ AdEq.CMA_DD.CMA_mu*ones(1,AdEq.CMA_DD.hardSwitch),...
    AdEq.CMA_DD.muInit*AdEq.CMA_DD.muFactor.^(min((AdEq.CMA_DD.muShifts+1),...      % vector of step size
    ceil((AdEq.CMA_DD.hardSwitch+1:Nsymb)/AdEq.CMA_DD.muInt))-1)];  
AdEq.CMA_DD.startmeas = AdEq.CMA_DD.muShifts*AdEq.CMA_DD.muInt+AdEq.CMA_DD.hardSwitch;

startmeas = max([Nsymb - 3000+1 , AdEq.LMS.startmeas, AdEq.CMA_DD.startmeas]);
if (startmeas ~= Nsymb - 3000+1)
    display('***** Warning: Wrong startmeas, better to use constant number of symbol for BER measurements');
    keyboard;
end


%% ======== Carrier Phase Recovery ========
Cpr.type        = 'DD';                                                    % carrier recovery type. 'dd' = decision-directed; 'nd' = non-data-aided
Cpr.Filter.Type = 'FIR';                                      
Cpr.Filter.Lmax = 40;
Cpr.Filter.f    = 0.95;


%% ============ Phase Tracker =============
Pt.mu = .005;

% Training mode vs. decision-directed mode ----------------------------
Pt.startdd  = AdEq.LMS.startdd;                                                          % symbol index at which to stop using training symbols and change to decision-directed mode
Pt.decwtvec = (1:Nsymb)<Pt.startdd;                                         % vector of decision weights (1 when use training symbols, 0 when use decision symbols)


%% ============ Launched Power ============
PlaunchdBm = -1;                        
PindBm = PlaunchdBm + DROP.ILdB;      
Pin = 1e-3*10^(PindBm/10);

%% ========= Performance Metrics =========
sizeAtt = size(Att.att);
Prec    = zeros(sizeAtt);                                                   % Signal power (W)
Srec    = zeros(sizeAtt);                                                   % Noise PSD (W/Hz)
OSNR    = zeros(sizeAtt);                                                   % OSNR per carrier
MSE_AWGN_expected = zeros(M,sizeAtt(2));                                    % Expected MSE at equalizer output
MSE_meas          = zeros(M,sizeAtt(2));                                    % Measured MSE at equalizer output
% Averaged over I and Q in two polarizations ------------------------------
Qmeas             = zeros(length(K),sizeAtt(2));                            % Measured Q factor 
BERmeas           = zeros(length(K),sizeAtt(2));                            % Measured BER 
BER_from_Qmeas    = zeros(length(K),sizeAtt(2));                            % Expected BER 

%% ============== Plots ===============
Plots = containers.Map();                                                   % List of figures 
Plots('CarrierSeparation')    = 0;
Plots('ModElectrodeLoss')     = 0;
Plots('CarrierPhaseNoise')    = 0;
Plots('DiffGroupDelay')       = 0;
Plots('SuperChannelSpectrum') = 0;
Plots('Equalizer')            = 0;
Plots('EqualizerTaps')        = 0;
Plots('PhaseTracker')         = 0;
Plots('PerformanceMetrics')   = 0;
Plots('PerformancePlots')     = 0; 
[~,ind]    = sort(omega);
[~,inddsp] = sort(omegadsp(abs(omegadsp)<=(pi*oversamp*Rs)));

%% =============== Figures ===============       
% Carrier Separation ------------------------------------------------------
if (Plots('CarrierSeparation'))
    figure(1)
    hold on
    plot(omegadsp(inddsp)/(2*pi*1e9),abs(Hcarsep(inddsp)),'b-','LineWidth',1.5) 
    plot(omegadsp(inddsp)/(2*pi*1e9),abs(HT2nd(inddsp)),'b--',...
        omegadsp(inddsp)/(2*pi*1e9),abs(HT4nd(inddsp)),'b-.',...
        omegadsp(inddsp)/(2*pi*1e9),abs(HSINC(inddsp)),'r'); 
    hold off
    axis tight
    xlabel('Frequency (GHz)'); ylabel('Magnitude Response')
    title('Carrier separation digital filters')
    legend('Cascade of two \itT\rm/2 and one \itT\rm/4','\itT\rm/2 delay-and-add', '\itT\rm/4 delay-and-add','Sinc function')
end
% Electorde Loss ----------------------------------------------------------
if (Plots('ModElectrodeLoss'))
    figure(2)
    plot(omega(ind)/2/pi*1e-9,10*log10(abs(Mod.Hel(ind))),'.');
    axis([-Rs/1e9*(M+1)/2 Rs/1e9*(M+1)/2 -5 0])
    xlabel('Frequency (GHz)'); ylabel('Magnitude Response')
    title('Electrode Loss frequency response')
end
% Differential Group Delay ------------------------------------------------
if (Plots('DiffGroupDelay'))
    figure(3)
    tau = zeros(1,length(omega)-1);
    for m = 1:length(omega)-1;
        tau(m) = 2/(omega(2)-omega(1))*sqrt(det(Fiber.M(:,:,m+1)-Fiber.M(:,:,m)));
    end
    plot(omega(1:end-1)/2/pi/1e9, abs(tau)/1e-12, 'b.');
    axis([-Rs/1e9*(M+1) Rs/1e9*(M+1) 0 20])
    xlabel('Frequency \omega/2\pi (GHz)'); ylabel('|\Delta\tau(\omega)| (ps)');
    title('Differential group delay vs. \omega');
end
%% ============= Modulation ==============
b = prbs(PRBS.N);%,PRBS.tabs);                                                 % Binary bit sequence (0,1)
Shifts = floor(Nsymb*(0:2*log2(Const.M)*M*(1+2*nSC)-1)/(2*log2(Const.M)*M*(1+2*nSC)));
if (mod(Ncar,M)) 
    error('Number of Carriers in Superchannel Should be Multiple of M');
end

c = ['b', 'y', 'r', 'g', 'm', 'c', 'b', 'y', 'r', 'g', 'm', 'c'];
display('*** Detecting ...') 
for k=K
    superChannel =  k;                                                      % Indices of M carriers being detected    
    subCarriers  =  (1:M)+(superChannel-1)*M;
    neighborSC   =  (superChannel-1-nSC)+1:(superChannel+nSC);
    validSuperCh = (neighborSC>=1)&(neighborSC<=Ncar/M);      
    fLO = fcars(subCarriers(1)) + (M-1)*Rs/2;                               % Frequency of local oscillator w.r.t. (Hz)   
    Ein = zeros(2,length(t));
    xdet = zeros(2*M,Nsymb);
    for sCh = 1:1+2*nSC
        if (validSuperCh(sCh))
            ind = 2*log2(Const.M)*M*(sCh-1)+1:2*log2(Const.M)*M*sCh;
            [Vin_SC, Xout] = QAM_SC_Tx(b,Shifts(ind),Tx,Rs,M);
            if (neighborSC(sCh)==superChannel) 
                Tx_PN = PN.Theta(sCh,:);
                xdet = Xout;
            end
            
            %filter drive waveforms for modulators
            Vix = real(ifft(fft(Vin_SC(1,:)).*Mod.Hel));
            Vqx = real(ifft(fft(Vin_SC(2,:)).*Mod.Hel)); 
            Viy = real(ifft(fft(Vin_SC(3,:)).*Mod.Hel)); 
            Vqy = real(ifft(fft(Vin_SC(4,:)).*Mod.Hel));
            % modulate signal fields (each has unit average power)
            E1tx   = sin(pi*Vix/2/Tx.Vpi) + 1i*sin(pi*Vqx/2/Tx.Vpi);        % x polarization
            E2tx   = sin(pi*Viy/2/Tx.Vpi) + 1i*sin(pi*Vqy/2/Tx.Vpi);        % y polarization
            
            EmpScale = 2/sqrt(mean(abs(E1tx).^2+abs(E2tx).^2));

            Ein_SC = EmpScale*sqrt(M*Pin)/2* ...                            % polarization multiplexed signal
               [exp(1i*PN.Theta(sCh,:)).*E1tx; exp(1i*PN.Theta(sCh,:)).*E2tx]; 
            
            if (mean(abs(Ein_SC(1,:)).^2+abs(Ein_SC(2,:)).^2)-(M*Pin))>.01*M*Pin
                display('***** Error: Wrong Input Power');
                keyboard;
            end
            
            delta_f = ((pattern == 1 && sCh == (nSC+1)) + (pattern == 2 && sCh <= nSC)- ...
                (pattern == 2 && sCh >= nSC+2))*df;
            if (FreqError == 0 && delta_f~=0)
               display('**** Error: Non-zero frequency')
               keyboard;
            end
            fSC    = fcars(neighborSC(sCh)*M-M+1) + (M-1)*Rs/2;             % Frequency of nth carrier(Hz)
            fIFn   = fSC - fLO + delta_f;                                   % Intermediate frequency of nth carrier in group (frequency w.r.t. local oscillator) (Hz)          

            Ein    = Ein + freqshift(Ein_SC,t,fIFn);
            %2/sqrt(mean(abs(E1tx).^2+abs(E2tx).^2))
            if Plots('SuperChannelSpectrum') 
                figure(4)
                hold on
                temp = freqshift(Ein_SC,t,fIFn);
                temp = smooth(abs(fft(temp(1,:))),16);
                plot(omega/(2*pi*1e9),10*log10(temp),'-',...
                    'MarkerEdgeColor',c(sCh),'MarkerSize',4);
                xlabel(' Frequency (GHz)'); ylabel('Magnitude Response')
                title('Super Channel Spectrum')
                axis([-60 60 5 20])
                %axis([-5*Rs/1e9*M 5*Rs/1e9*M 5 20])
                hold off
            end
        end
    end

%% ========= Network Propagation ========== 
     Eout = Ein;                                                            % Initial electric field waveform
     Pout = Pin;                                                            % Initial value of nominal output signal power (W)
     Sout = 0;                                                              % Initial value of nominal output noise PSD (W/Hz)

    % Add -------------------------
    % First Amplifier -------------
    G1nom = 10^(Amp1.GdB/10);                                               % Nominal power transmittance of second amplifier (W/W)
    S1nom = Amp1.nsp*(G1nom-1)*h*nu;                                        % Nominal PSD of ASE from second amplifier (W/Hz)
    G1 = G1nom*ones(size(omega));                                           % Frequency-dependent power transmittance of second amplifier (W/W)
    S1 = Amp1.nsp*(G1-1)*h*nu;                                              % Frequency-dependent PSD of ASE from second amplifier (W/Hz)
    % Fiber ----------------------- 
    % Second amplifier ------------
    G2nom = 10^(Amp2.GdB/10);                                               % Nominal power transmittance of second amplifier (W/W)
    S2nom = Amp2.nsp*(G2nom-1)*h*nu;                                        % Nominal PSD of ASE from second amplifier (W/Hz)
    G2 = G2nom*ones(size(omega));                                           % Frequency-dependent power transmittance of second amplifier (W/W)
    S2 = Amp2.nsp*(G2-1)*h*nu;                                              % Frequency-dependent PSD of ASE from second amplifier (W/Hz)

    % Drop ------------------------
    
    for span=1:Fiber.Nspan
        Pout = Pout*Add.T*G1nom*Fiber.Tf*G2nom*Drop.T;                          % nominal signal power at output of segment
        Sout = Sout*Add.T*G1nom*Fiber.Tf*G2nom*Drop.T + S1nom*Fiber.Tf*G2nom*Drop.T + S2nom*Drop.T;

        % propagate signals through segment
        Eadd = scprop(Eout,Add.H);                                              % electric field at output of add module
        Eamp1 = scamp(Eadd,deltat,G1,S1);                                       % electric field at output of first amplifier
        Ef = scprop(Eamp1,Fiber.H);                                             % electric field at output of fiber
        Ef = vprop(Ef,Fiber.M);
        Eamp2 = scamp(Ef,deltat,G2,S2);                                         % electric field at output of second amplifier
        Eout = scprop(Eamp2,Drop.H);                                            % electric field at output of drop module = output of segment
    end
    
    display(subCarriers)
    
%% ========== Nosie Loading ============== 
    for n=1:length(Att.att) 
        att=Att.att(n);
        % attenuate and amplify signal ------------------------------------ 
        Gnom = 10^abs(att/10);                                              % Gain of noise loading amplifier (W/W)
        G = Gnom*ones(size(omega));
        Snom = Amp1.nsp*(Gnom-1)*h*nu;                                      % PSD of ASE added by noise loading amplifier (W/Hz)
        S = Snom*ones(size(omega)); 
        Prec(n) = Pout;                                                  % Received signal power (W)
        Srec(n) = Sout + Snom;                                             % Received PSD of ASE (W/Hz)
        Eat = scprop(Eout,1./sqrt(G));                                      % Attentuate
        Erec = scamp(Eat,deltat,G,S);                                       % Amplify
        % OSNR ------------------------------------------------------------
        osnrnoisebw = 12.5e9;                                               % Optical filter BW for measurement of optical signal-to-noise ratio
        OSNR(n) = Prec(n)/(2*Srec(n)*osnrnoisebw);                          % Received OSNR

        display([' @ OSNR = ' num2str(10*log10(OSNR(n)))]);
%% ============== Receiver ================
        % Coherent detection ----------------------------------------------
        [Y agcgain] = PDM_QAM_Rx(Erec,xdet,Rx,t,omega);
        % Sampling and quantization ---------------------------------------
        [Ys tsamp]=Sampler(t, Y, Nstep, Nsymb, oversamp);
        Ysq = Quant(Ys,Quantz);
        % CD compensation ------------------------------------------------- 
        Ycdcomp = scprop(Ysq,Hcdcomp); 
%% ======== Digital signal processing =========
        for l = 1:M                                                         % Index if subcarrier
            m = subCarriers(l);                                             % index of carrier being detected
            % Shift carrier to baseband -----------------------------------
            % # NOT NEEDED FOR SINGLE CARRIER 
            Yshift = freqshift(Ycdcomp,tdsp,(M-2*l+1)*Rs/2+df_LO);
            % Carrier separation ------------------------------------------
            % # NOT NEEDED FOR SINGLE CARRIER    
            Ysep = scprop(Yshift,Hcarsep);                                  % signal after carrier separation
            % Adaptive equalization ---------------------------------------
            PLO = 1e-3*10^(Rx.PLOdBm/10); 
            if strcmp(Rx.PolSplit.sig,'PBS')&&strcmp(Rx.PolSplit.LO,'PBS');
                MSE_AWGN_expected(l,n)=2*Srec(n)*(PLO/2)*Rs*agcgain^2;      % rough estimate of expected MSE at equalizer output (unitless)
            else
                MSE_AWGN_expected(l,n)=2*Srec(n)*(PLO/2)*Rs*agcgain^2/2;    % rough estimate of expected MSE at equalizer output (unitless)
            end
            x = xdet(2*l-1:2*l,:);  
            OSNR_int = int32(10*log10(OSNR(n)));
            if (DumpData == 1 && OSNR_int==SNR)
                outfile      = sprintf('M%d_Car%d_OSNR%d_nSC%d.mat',M,m,OSNR_int,nSC);
                save(outfile,'x','tsamp','Ysep','Nsymb','oversamp','PN');
                keyboard;
            end
            
            switch Equalization
                case 'CMA'       
                    [xequ, eps_cma, WT] = CMA(tsamp, Ysep, Nsymb, oversamp, AdEq.CMA);
                    if (PhaseNoise == 1)
                        [xcpr, ~, thetahat, Delta] = CPR(xequ,Cpr,PN,OSNR(n));  
                        [xD, xhat, eps, theta_pt] = PT(x,xcpr, Nsymb,Pt);
                        measind = (startmeas:Nsymb)-Delta;                  % indices of symbols during which to measure 
                    else
                        [xD, xhat, eps, theta_pt] = PT(x,xequ, Nsymb,Pt);
                        measind = (startmeas:Nsymb);                        % indices of symbols during which to measure 
                    end   
                    if (OSNR_int==SNR && k==1)
                        plot_CMA;
                        plot_EquTaps
                        plot_CPR; 
                        plot_PT;
                    end
                case 'MMA'
                     [xhat, eps, WT] = MMA(tsamp, Ysep, Nsymb, oversamp, AdEq.MMA); 
                     xD = -sign(real(xhat))-1i*sign(imag(xhat));
                     measind = (startmeas:Nsymb);
                     %plot_MMA;
                case 'LMS'
                    if (PhaseNoise == 1)
                        display('***** Error: LMS equalizer alone does not work in the presence of phase noise');
                        keyboard;
                    end
                    [xhat, xD, eps, ~, WTe, WTo] = LMS(x, tsamp, Ysep, Nsymb, oversamp, AdEq.LMS);
                    measind = startmeas:Nsymb;
                    if (OSNR_int==SNR && k==1)
                        plot_LMS
                    end

                case 'CMA-DD'
                    if (PhaseNoise == 1)
                        [xhat, xD, eps, WT, thetahat, Delta] = CMA_DD_CPR(tsamp, Ysep, Nsymb, oversamp, AdEq.CMA_DD, Cpr , PN, SNR);
                        measind = startmeas:Nsymb-Delta;
                        xcpr= xhat;
                    else
                        [xhat, xD, eps, WT] = CMA_DD(tsamp, Ysep, Nsymb, oversamp, AdEq.CMA_DD);
                        measind = startmeas:Nsymb;
                    end
                    if (OSNR_int==SNR && k==6)
                        plot_CMA_DD
                        plot_CPR
                        %plot_EquTaps
                    end
                    [xD, xhat, eps, theta_pt] = PT(x,xhat, Nsymb,Pt);
                    if (OSNR_int==SNR && k==1)
                        plot_PT;
                    end
                otherwise
                    display('***** Error: Unrecognized Equalizer')
            end;
            
            % Performance metrics -----------------------------------------
            % Compute equalizer error, Q and BER
            numsymbmeas = length(measind);                                  % number of symbols measured
            x1imeas = real(x(1,measind)); x1qmeas = imag(x(1,measind));     % correct symbols over indices measured
            x2imeas = real(x(2,measind)); x2qmeas = imag(x(2,measind));
            % estimated symbols over indices measured
            x1ihatmeas = real(xhat(1,measind)); x1qhatmeas = imag(xhat(1,measind)); 
            x2ihatmeas = real(xhat(2,measind)); x2qhatmeas = imag(xhat(2,measind));          
            % decision symbols over indices measured                       
            x1iDmeas = real(xD(1,measind)); x1qDmeas = imag(xD(1,measind));       
            x2iDmeas = real(xD(2,measind)); x2qDmeas = imag(xD(2,measind));
            % true error for indices measured
            eps1imeas = real(eps(1,measind)); eps1qmeas = imag(eps(1,measind));     
            eps2imeas = real(eps(2,measind)); eps2qmeas = imag(eps(2,measind));
            sqerr = abs(eps(1,:)).^2 + abs(eps(2,:)).^2;
            MSE_meas(l,n) = mean(sqerr(measind));
            % Q factor (based on true error, starting at 'startmeas')
            Q1imeas = (mean(x1ihatmeas(x1imeas>0))-mean(x1ihatmeas(x1imeas<0)))/(std(eps1imeas(x1imeas<0))+std(eps1imeas(x1imeas>0)));
            Q1qmeas = (mean(x1qhatmeas(x1qmeas>0))-mean(x1qhatmeas(x1qmeas<0)))/(std(eps1qmeas(x1qmeas<0))+std(eps1qmeas(x1qmeas>0)));
            Q2imeas = (mean(x2ihatmeas(x2imeas>0))-mean(x2ihatmeas(x2imeas<0)))/(std(eps2imeas(x2imeas<0))+std(eps2imeas(x2imeas>0)));
            Q2qmeas = (mean(x2qhatmeas(x2qmeas>0))-mean(x2qhatmeas(x2qmeas<0)))/(std(eps2qmeas(x2qmeas<0))+std(eps2qmeas(x2qmeas>0)));           
            Qmeas(m,n) = (Q1imeas+Q1qmeas+Q2imeas+Q2qmeas)/4;
            % Bit-error ratio, average over all four tributaries, starting at 'startmeas', computed based on measured Q factors
            BER_from_Qmeas(m,n) = (0.5*erfc(Q1imeas/sqrt(2))+0.5*erfc(Q1qmeas/sqrt(2))+...
                0.5*erfc(Q2imeas/sqrt(2))+0.5*erfc(Q2qmeas/sqrt(2)))/4;
            % measured
            num_bit_err_meas = sum(x1iDmeas~=x1imeas)+sum(x1qDmeas~=x1qmeas)+sum(x2iDmeas~=x2imeas)+sum(x2qDmeas~=x2qmeas);
            BERmeas(m,n) = num_bit_err_meas/4/numsymbmeas;
                        % Plot equalizer output         
                        
%             keyboard; 
%             figure
%             for k=measind(1):measind(100)
%                 plot(x(2,k),'ok','MarkerEdgeColor','b','MarkerSize',7)
%                 axis([-2 2 -2 2])
%                 hold on
%                 plot(xequ(2,k),'sr','MarkerEdgeColor','r','MarkerSize',10)
%                 %plot(xhat(2,k),'ok','MarkerEdgeColor','g','MarkerSize',10)
%                 %plot(xD(2,k),'*k','MarkerEdgeColor','b','MarkerSize',5)
%                 pause(1)
%                 hold off
%             end  
% 
%             
            % List performance metrics
            if Plots('PerformanceMetrics')
                figure
                axis off                
                text(0,1,'Receiver Filter');
                text(0,0.95,['     Type = ' Rx.LPF.Type ', Order = ' num2str(Rx.LPF.Ord) ', Bandwidth = ', num2str(Rx.LPF.BW/Rs), 'Rs']);
                text(0,0.8,['OSNR (per carrier in ' num2str(osnrnoisebw/1e9,3) ' GHz) = ' num2str(10*log10(OSNR(n)),3) ' dB']);
                text(0,0.7,['Carrier ' num2str(m) ', ' num2str(oversamp) ' samples/symbol']);
                text(0,0.6,['10 log_{10}(\itQ_{meas}\rm^{2}) = ' num2str(10*log10(Qmeas(m,n)^2),3) ' dB']);
                text(0,0.5,['BER calc. from \itQ_{meas}\rm = ' num2str(BER_from_Qmeas(m,n),3)]);
                text(0,0.4,['BER_{meas} = ' num2str(BERmeas(m,n),3) ', num. bit err. meas. = ' num2str(num_bit_err_meas,3)]);
            end
        end
    end
end

%keyboard;   
out.M        = M;
out.Nsymb    = Nsymb;
out.AdEq     = AdEq;
out.OSNR     = OSNR;
out.Qmeas    = Qmeas;
out.BERmeas  = BERmeas;
out.BER_from_Qmeas = BER_from_Qmeas;
if (PMD==1)
    outfile = sprintf('results/NGI_OFDM/PMD/noPN/%s/M%d/nSC%d',Equalization,M,nSC);
else
    if (PhaseNoise == 1)
        outfile = sprintf('results/NGI_OFDM/noPMD/PN%d/%s/M%d/nSC%d',PN.d_nu/1e3,Equalization,M,nSC);
    else
        outfile = sprintf('results/NGI_OFDM/noPMD/noPN/%s/M%d/nSC%d/run%d',Equalization,M,nSC,run);
    end
end
if (FreqError == 1)
    outfile = sprintf('results/NGI_OFDM/noPMD/noPN/TxFreqErr/%s/M%d/nSC%d/ratio%d',Equalization,M,nSC, error_ratio);
end

if(Save==1)
    if (~exist(outfile,'dir'))
        mkdir(outfile);
        mkdir([outfile '/fig']);
    end
    save([outfile '/result'],'OSNR','Qmeas','BERmeas','BER_from_Qmeas');
    display('Saved')
    
    H = get(0,'children');
    for n = 1: length(H)
        saveas(H(n),[outfile '/fig/' num2str(H(n)) '.fig']);
        close(H(n));
    end
end


% COMPUTE AND PLOT PERFORMANCE METRICS FOR ALL CARRIERS
% Compute average performance metrics
if Plots('PerformancePlots') 
    %geommean_of_Qmeas = exp(mean(log(Qmeas),1));                                % Geometric mean over all carriers of Qmeas
    %arithmean_of_BER_from_Qmeas = mean(BER_from_Qmeas,1);                       % Arithmetic over all carriers of BER computed from Qmeas
    %BER_from_geommean_of_Qmeas = 0.5*erfc(geommean_of_Qmeas/sqrt(2));           % BER computed from geometric mean over all carriers of Qmeas
	
    %Plot Q factor
% 	figure
% 	plot(10*log10(OSNR),10*log10(Qmeas.^2)','g--',10*log10(OSNR),10*log10(geommean_of_Qmeas.^2)','k-');
% 	g=get(gca,'children'); set(g,'linewidth',1.5);
% 	xlabel(['OSNR per carrier in ' num2str(osnrnoisebw/1e9,3) '-GHz BW (dB)']); ylabel('10 log_{10}(\itQ\rm_{meas}\rm^{2}) (dB)');
% 	
% 	figure
%     semilogy(10*log10(OSNR),BER_from_Qmeas','g-',10*log10(OSNR),BER_from_geommean_of_Qmeas','k-');
% 	g=get(gca,'children'); set(g,'linewidth',1.5);
% 	xlabel(['OSNR per carrier in ' num2str(osnrnoisebw/1e9,3) '-GHz BW (dB)']); ylabel('BER');
%     %legend('Measured','Expected from Q_{meas}');
%     %axis([10 20 1e-7 1e-2])
	figure
    semilogy(10*log10(OSNR),BER_from_Qmeas','g-');
	g=get(gca,'children'); set(g,'linewidth',1.5);
	xlabel(['OSNR per carrier in ' num2str(osnrnoisebw/1e9,3) '-GHz BW (dB)']); ylabel('BER');
    %legend('Measured','Expected from Q_{meas}');
    %axis([10 20 1e-7 1e-2])

%     
%     figure
%     semilogy(10*log10(OSNR),BER_from_Qmeas','g-',10*log10(OSNR),BERmeas(1,:),'o')
%     g=get(gca,'children'); set(g,'linewidth',1.5);
% 	xlabel(['OSNR per carrier in ' num2str(osnrnoisebw/1e9,3) '-GHz BW (dB)']); ylabel('BER');
%     
    figure
    axis off                
    text(0,1, 'System:')
    text(.1,.9,['CD: 1     PMD: ' num2str(PMD) '     PN: ' num2str(PhaseNoise)]) 
    text(.1,.8,['Equalization: ' Equalization])
    text(.1,.7,['Carrier Phase Recovery: ' Cpr.type])
    text(0,.6,'Rx Filter');
    text(.1,0.5,[' Type : ' Rx.LPF.Type ', Order : ' num2str(Rx.LPF.Ord) ', Bandwidth : ', num2str(Rx.LPF.BW/Rs), 'Rs']);
    text(0,.4,['Adaptive Equalizer: ' Equalization]);
    if (strcmp(Equalization,'LMS'))
        text(.1,.3,['L_{Filt} : ' num2str(AdEq.LMS.LFilt)])
        text(.1,.2,['\mu_{initial} : ' num2str(AdEq.LMS.muInit)])
        text(.1,.1,['\mu_{factor} : ' num2str(AdEq.LMS.muFactor)])
        text(.1,0,['\mu_{interval : }' num2str(AdEq.LMS.muInt)])
        text(.1,-.1,['Number of shifts : ' num2str(AdEq.LMS.muShifts)])
    elseif (strcmp(Equalization,'CMA'))
        text(.1,.3,['L_{Filt} : ' num2str(AdEq.CMA.LFilt)])
        text(.1,.2,['\mu : ' num2str(AdEq.CMA.mu)])
    end
if (Save == 1)
     H = get(0,'children');
     for n = 1: length(H)
         saveas(H(n),[outfile '/' num2str(H(n)) '.fig']);
         close(H(n));
     end
end
end


           

