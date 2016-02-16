% Ngi_CoOFDM.m     Milad Sharif
% No-Gaurd Interval Co-OFDM with M subcarriers per sampling and one
% modulator for group of M subcarriers.
% PMD, PDL, higher-order modulation, laser phase noise, fiber nonlinearity
% are not included 
% thermal noise, LO shot noise are neglected.
clear, clc, close all

addpath plots_scripts/
addpath DSP/
addpath f/
addpath ../f/
addpath ../apd/
addpath ../soa/

%% ============== Simulation parameters =====================
sim.Nsymb = 2^15; % Number of symbols in montecarlo simulation
sim.ros = 5/4;
sim.Mct = 8*sim.ros;    % Oversampling ratio to simulate continuous time 
sim.BERtarget = 1.8e-4; 
sim.Ndiscard = 64; % number of symbols to be discarded from the begining and end of the sequence (should be larger than Ntrain, Ntaps, Nfft, etc
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation
sim.Rb = 100e9; % Bit rate
sim.M    = 4;                                                              % QAM order
sim.pulse = 'NRZ';                                                         % Transmitter pulse shape ('NRZ')
sim.Npol = 2;
sim.Rs = sim.Rb/(sim.Npol*log2(sim.M));

% Simulation
sim.RIN = true; 
sim.PMD = ~true;
sim.phase_noise = true;
sim.pre_amplification = false;
sim.quantiz = false;

%% Time and frequency
sim.fs = sim.Rs*sim.Mct;  % sampling frequency in 'continuous-time'

dt = 1/sim.fs;
t = 0:dt:(sim.N-1)*dt;
df = 1/(dt*sim.N);
f = -sim.fs/2:df:sim.fs/2-df;

sim.t = t;
sim.f = f;

%% =========== Mod. ==================
omega = 2*pi*f;
c = 3e8;
Mod.Vpi = 4.5;                                                             % Modulator switching voltage (V)
ratio      = 1; %.96;                                                            
n_r        = 2.15;                                                          
n_m        = n_r*ratio; 
Mod.d_12   = (n_m -n_r)/c;                                                  % Walkoff parameter between microwave and optical fields
Mod.a      = .005*100*sqrt(abs(omega)/2/pi*1e-9);                            % Microwave attenuation coefficient (electrode loss)
Mod.L      = .05;                                                           % Modulator length
Mod.Hel    = (1-exp(Mod.L*(-Mod.a + 1j*Mod.d_12*omega)))./...               % Freq. response of the optical modulator limited by elec. loss and velocity mismatch
                (Mod.a-1j*Mod.d_12*omega)/Mod.L;
Mod.Hel(isnan(Mod.Hel)) = 1;                                               % Mod.Hel(f=0) is NaN            
%% =========== Transmitter ========== 
Tx.Laser.PdBm = 0;                                                         % output power in dBm
Tx.Laser.lamb = 1250e-9;                                                   % wavelength (m) and reference frequency
Tx.Laser.RIN = -150;                                                       % relative intensity noise dB/Hz
Tx.Laser.linewidth = 200e3;                                                % linewidth (Hz)
Tx.Laser.freqshit = 0;                                                     % frequency shift with respect to reference frequency

Tx.Mod = Mod;                                                              % optical modulator

Tx.filt = design_filter('bessel', 5, 1.2*sim.Rs/(sim.fs/2));                    % design_filter(type, order, normalized cutoff frequency)
Tx.Delx  = 0;                                                               % Delay of x pol. in pol. mux. (symbol intervals)
Tx.Dely  = 0;                                                               % Delay of y pol. in pol. mux. (symbol intervals)

%% ====== Fiber ========
%  Constructor: fiber(L, att(lamb) (optional), D(lamb) (optional)) 
% L : fiber length (m)
% att(lamb) : function handle of attenuation (att) at wavelength (lamb),
% deafault is att(lamb) = 0 dB/km
% D(lamb) : function handle of dispersion (D) at wavelength (lamb) in ps/(kmnm),
% default is D(lamb) = SSMF with lamb0 @ 1310 ps/(kmnm)
Fiber = fiber(10e3);
Fiber.PMD = sim.PMD;                                                       % whether to similate PMD
Fiber.meanDGDps = 4.5;                                                     % Mean DGD (ps)
Fiber.PMD_section_length = 1e3;                                            % Controls number of sections to simulate PMD (m)

%% =========== Amplifier ==============
% Constructor: soa(GaindB, NF, lamb, maxGain (optional))
% GaindB : Gain in dB
% NF : Noise figure in dB
% lamb : wavelength in m
% maxGain = maximum amplifier gain in dB, default is Inf
Amp = soa(20, 7, Tx.Laser.lamb);
% Note: class soa can be used for any amplifier, since it essentially 
% characterizes the amplifier in terms of gain and noise figure only

%% =============== Receiver ===============
% Local oscillator --------------------------------------------------------
Rx.LO = Tx.Laser;                                                          % Copy parameters from TX laser
RX.LO.PdBm = 13;                                                           % Total local oscillator power (dBm)
Rx.LO.freqshift = 1e-2*1e9;                                                    % Frequency shift with respect to transmitter laser in Hz
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
% Constructor: pin(R, Id, BW (optional))
% R : responsivity in A/W
% Id : dark current in A
% BW : bandwidth, default = Inf. Frequency response is a
% first-order filter with bandwidth BW.
Rx.PD = pin(1, 10e-9);

Rx.N0 = (30e-12)^2;

%% ============= ADC ==============
Rx.ADC.ENOB = 5;                                                           % effective number of bits
Rx.ADC.rclip = 0;                                                          % clipping ratio: clipped intervals: [xmin, xmin + xamp*rclip) and (xmax - xamp*rclip, xmax]
Rx.ADC.ros = sim.ros;                                                      % oversampling ratio with respect to symbol rate 
Rx.ADC.fs = Rx.ADC.ros*sim.Rs;                                             % ADC sampling rate
Rx.ADC.filt = design_filter('cheby1', 5, 0.5*Rx.ADC.fs/(sim.fs/2));        % design_filter(type, order, normalized cutoff frequency)

%% ========== Digital Filters ===== ======
% CD compensation ---------------------------------------------------------
% wintype = 'rectwin';                                                        % Window type. Options include: bartlett, chebwin, hamming, hanning, rectwin (see MATLAB "help window"). 
% eval(['win = window(@' wintype ', length(omegadsp));']);
% Hcdcomp = conj(Fiber.Hdisp(sim.f, Tx.lamb)).*win.';                        % Freq resp of CD compensating digital filter

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
    ceil((1:sim.Nsymb)/AdEq.LMS.muInt))-1);     
% Training mode vs. decision-directed mode ----------------------------
AdEq.LMS.startdd = round(AdEq.LMS.muShifts/2)*AdEq.LMS.muInt;               % symbol index at which to stop using training symbols and change to decision-directed mode
AdEq.LMS.decwtvec = (1:sim.Nsymb)<AdEq.LMS.startdd;                             % vector of decision weights (1 when use training symbols, 0 when use decision symbols)
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
    ceil((AdEq.CMA_DD.hardSwitch+1:sim.Nsymb)/AdEq.CMA_DD.muInt))-1)];  
AdEq.CMA_DD.startmeas = AdEq.CMA_DD.muShifts*AdEq.CMA_DD.muInt+AdEq.CMA_DD.hardSwitch;

startmeas = max([sim.Nsymb - 3000+1 , AdEq.LMS.startmeas, AdEq.CMA_DD.startmeas]);
if (startmeas ~= sim.Nsymb - 3000+1)
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
Pt.decwtvec = (1:sim.Nsymb)<Pt.startdd;                                         % vector of decision weights (1 when use training symbols, 0 when use decision symbols)


%% ============ Launched Power ============
PlaunchdBm = -1;                        
PindBm = PlaunchdBm;      
Pin = 1e-3*10^(PindBm/10);

%% ========= Performance Metrics =========
% sizeAtt = size(Att.att);
% Prec    = zeros(sizeAtt);                                                   % Signal power (W)
% Srec    = zeros(sizeAtt);                                                   % Noise PSD (W/Hz)
% OSNR    = zeros(sizeAtt);                                                   % OSNR per carrier
% MSE_AWGN_expected = zeros(M,sizeAtt(2));                                    % Expected MSE at equalizer output
% MSE_meas          = zeros(M,sizeAtt(2));                                    % Measured MSE at equalizer output
% % Averaged over I and Q in two polarizations ------------------------------
% Qmeas             = zeros(length(K),sizeAtt(2));                            % Measured Q factor 
% BERmeas           = zeros(length(K),sizeAtt(2));                            % Measured BER 
% BER_from_Qmeas    = zeros(length(K),sizeAtt(2));                            % Expected BER 

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

%% ============= Modulation ==============
% b = prbs(log2(sim.sim.Nsymb)); % PRBS
dataTX = randi([0 sim.M-1], [2, sim.Nsymb]); % symbol stream for each polarization

[Vin, symbolsTX] = QAM_SC_Tx(dataTX, Tx, sim); % generates QAM signal

Ein = cw_laser(Tx.Laser, sim); % Generates electric field with intensity and phase noise

Ein = eom(Ein, Vin, Tx.Mod); % modulate optical signal using eletro-optical modulator (EOM)

% ========= Propagation ========== 
Erec = Fiber.linear_propagation(Ein, sim.f, Tx.Laser.lamb);

% ============== Receiver ================
c = ['b', 'y', 'r', 'g', 'm', 'c', 'b', 'y', 'r', 'g', 'm', 'c'];
display('*** Detecting ...') 

Y = PDM_QAM_Rx(Erec, sim.M, Rx, sim);

% Digital to analog conversion: filters, sample, and quantize
[Yxi, varQ(1)] = adc(real(Y(1, :)), Rx.ADC, sim);
[Yxq, varQ(2)] = adc(imag(Y(1, :)), Rx.ADC, sim);
[Yyi, varQ(3)] = adc(real(Y(2, :)), Rx.ADC, sim);
[Yyq, varQ(4)] = adc(imag(Y(2, :)), Rx.ADC, sim);

Ys = [Yxi + 1j*Yxq; Yyi + 1j*Yyq];

scatterplot(Ys(1, 1e4:end-512))
scatterplot(Ys(2, 1e4:end-512))

% Ys = [Yxi + 1j*(Yxq - Yxi*(sum(Yxi.*Yxq)/sum(Yxi.^2)));...
%     Yyi + 1j*(Yyq - Yyi*(sum(Yyi.*Yyq)/sum(Yyi.^2)))];

eq.Ntrain = Inf;
eq.mu = 1e-2;
eq.Ntaps = 5;
eq.ros = sim.ros;
[Y, Wx, Wy, MSE] = adaptive_td_fse(eq, Ys, symbolsTX, sim.M);
% [Y, Wx, Wy, MSE] = fse_cma(eq, Ys);
% [Y, Wx, Wy, MSE] = fse_cma_phase_track(eq, Ys, symbolsTX, sim.M);
% [Y, Wx, Wy, MSE] = fse_cma_dd(eq, Ys, sim.M);

scatterplot(Y(1, 1.5e4:end-512))
scatterplot(Y(2, 1.5e4:end-512))

% ============== Receiver ================
% ======== Digital signal processing =========
%         Adaptive equalization ---------------------------------------
switch upper(Equalization)
    case 'CMA'       
        [xequ, eps_cma, WT] = CMA(tsamp, Ysep, sim.Nsymb, oversamp, AdEq.CMA);
        if (PhaseNoise == 1)
            [xcpr, ~, thetahat, Delta] = CPR(xequ,Cpr,PN,OSNR(n));  
            [xD, xhat, eps, theta_pt] = PT(x,xcpr, sim.Nsymb,Pt);
            measind = (startmeas:sim.Nsymb)-Delta;                  % indices of symbols during which to measure 
        else
            [xD, xhat, eps, theta_pt] = PT(x,xequ, sim.Nsymb,Pt);
            measind = (startmeas:sim.Nsymb);                        % indices of symbols during which to measure 
        end   
        if (OSNR_int==SNR && k==1)
            plot_CMA;
            plot_EquTaps
            plot_CPR; 
            plot_PT;
        end
    case 'MMA'
         [xhat, eps, WT] = MMA(tsamp, Ysep, sim.Nsymb, oversamp, AdEq.MMA); 
         xD = -sign(real(xhat))-1i*sign(imag(xhat));
         measind = (startmeas:sim.Nsymb);
         plot_MMA;
    case 'LMS'
        if (PhaseNoise == 1)
            display('***** Error: LMS equalizer alone does not work in the presence of phase noise');
            keyboard;
        end
        [xhat, xD, eps, ~, WTe, WTo] = LMS(x, tsamp, Ysep, sim.Nsymb, oversamp, AdEq.LMS);
        measind = startmeas:sim.Nsymb;
        if (OSNR_int==SNR && k==1)
            plot_LMS
        end

    case 'CMA-DD'
        if (PhaseNoise == 1)
            [xhat, xD, eps, WT, thetahat, Delta] = CMA_DD_CPR(tsamp, Ysep, sim.Nsymb, oversamp, AdEq.CMA_DD, Cpr , PN, SNR);
            measind = startmeas:sim.Nsymb-Delta;
            xcpr= xhat;
        else
            [xhat, xD, eps, WT] = CMA_DD(tsamp, Ysep, sim.Nsymb, oversamp, AdEq.CMA_DD);
            measind = startmeas:sim.Nsymb;
        end
        if (OSNR_int==SNR && k==6)
            plot_CMA_DD
            plot_CPR
            plot_EquTaps
        end
        [xD, xhat, eps, theta_pt] = PT(x,xhat, sim.Nsymb,Pt);
        if (OSNR_int==SNR && k==1)
            plot_PT;
        end
    otherwise
        display('***** Error: Unrecognized Equalizer')
end

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
measured
num_bit_err_meas = sum(x1iDmeas~=x1imeas)+sum(x1qDmeas~=x1qmeas)+sum(x2iDmeas~=x2imeas)+sum(x2qDmeas~=x2qmeas);
BERmeas(m,n) = num_bit_err_meas/4/numsymbmeas;
%             Plot equalizer output         

keyboard; 
figure
for k=measind(1):measind(100)
    plot(x(2,k),'ok','MarkerEdgeColor','b','MarkerSize',7)
    axis([-2 2 -2 2])
    hold on
    plot(xequ(2,k),'sr','MarkerEdgeColor','r','MarkerSize',10)
    %plot(xhat(2,k),'ok','MarkerEdgeColor','g','MarkerSize',10)
    %plot(xD(2,k),'*k','MarkerEdgeColor','b','MarkerSize',5)
    pause(1)
    hold off
end  


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


%keyboard;   
out.M        = M;
out.sim.Nsymb    = sim.Nsymb;
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


           

