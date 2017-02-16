%% Validate ofdm class

%% Evaluation of OFDM in IM-DD system, which may be amplified or not
clear

% profile on

addpath f/
addpath ../f/
addpath ../apd/

tic

%% Simulation parameters
sim.Rb = 112e9;    % bit rate in bits/sec
sim.Nsymb = 2^11; % Number of symbols in montecarlo simulation
sim.Mct = 5;      % Oversampling ratio to simulate continuous time. Must be integer multiple of sim.ros.txDSP and numerator of sim.ros.rxDSP
sim.BERtarget = 1e-3; 
sim.Ndiscard = 256; % number of symbols to be discarded from the begining and end of the sequence
sim.OFDM = 'ACO-OFDM'; % {'DC-OFDM', 'ACO-OFDM'}
sim.quantiz = false;

%% OFDM 
% OFDM constructor: ofdm(Nc, Nu, CS, Rb, power_allocation_type (optional, default = 'palloc')
% Nc : number of subcarriers
% Nu : number of used subcarriers
% CS : nominal constellation size
% Rb : bit rate (b/s)
% power_allocation_type : {'palloc', 'preemphasis'}
if strcmpi(sim.OFDM, 'ACO-OFDM')
    disp('-- ACO-OFDM simulation')
    ofdm = ofdm(256, 208, 64, sim.Rb, 'palloc'); 
    ofdm.aco_ofdm_config();
elseif strcmpi(sim.OFDM, 'DC-OFDM')
    disp('-- DC-OFDM simulation')
    ofdm = ofdm(256, 208, 16, sim.Rb, 'palloc'); 
else
    error('sim.OFDM must be either DC-OFDM or ACO-OFDM')
end
ofdm.set_cyclic_prefix(5, 5); % set cyclic prefix length. Should be consistent with channel memory length

%% Time and frequency
sim.N = sim.Mct*(ofdm.Nc + ofdm.Ncp)*sim.Nsymb;           % total number of points simulated in continuous time
sim.fs = ofdm.fs*sim.Mct;
[sim.f, sim.t] = freq_time(sim.N, sim.fs);

%% ==========================  Transmitter ================================
%% DAC
DAC.fs = ofdm.fs; % DAC sampling rate
DAC.ros = 1; % oversampling ratio of transmitter DSP
DAC.resolution = 5; % DAC effective resolution in bits
DAC.filt = design_filter('butter', 5, 0.5*ofdm.fs/(sim.fs/2)); % DAC analog frequency response

%% ============================== Channel =================================
F = ClassFilter('butter', 5, 40e9/(sim.fs/2), sim.fs);
% F = ClassFilter('bessel', 5, 25e9/(ofdm.fs/2), ofdm.fs);
N0 = (30e-12).^2;

%% Receiver DSP
%% ADC
ADC.ros = 1;
ADC.fs = ofdm.fs;
ADC.filt = design_filter('butter', 5, 0.5*ofdm.fs/(sim.fs/2)); % Antialiasing filter
ADC.ENOB = 5; % effective number of bits. Quantization is only included if sim.quantiz = true and ENOB ~= Inf
% Rx.ADC.rclip = 0.05;

%% Equalizer
AdEq.mu = 1e-3;
AdEq.Ntrain = 512; % Number of frames used in training (if Inf all symbols are used)
% Note: training symbols are not used to compute the BER. Hence sim.Nsymb -
% Rx.AdEq.Ntrain must be large to obtain accurate BER estimate

%% Power allocation
noiseBW = ofdm.fs/2;
varNoise = 1/ofdm.Nc*N0*noiseBW*abs(ADC.filt.H(ofdm.fc/sim.fs)).^2;

Nhold = sim.Mct;
hZOH = 1/Nhold*ones(1, Nhold);
HZOH = freqz(hZOH, 1, ofdm.fc, sim.fs);
Hch = ADC.filt.H(ofdm.fc/sim.fs).*F.H(ofdm.fc).*DAC.filt.H(ofdm.fc/sim.fs).*HZOH;

ofdm.power_allocation(Hch, varNoise, sim.BERtarget, true);
% ofdm.power_allocation(ones(size(ofdm.fc)), varNoise, sim.BERtarget, true);

%% Generate OFDM signal
[xd, AdEq.trainSeq] = ofdm.signal(sim.Nsymb); 

%% DAC
xt = dac(xd, DAC, sim); % digital-to-analog conversion 

%% Channel
yt = F.filter(xt) + sqrt(N0*sim.fs/2)*randn(size(xt));

%% ADC
ADC.timeRefSignal = xt;
yk = adc(yt, ADC, sim);

% yk = F.filter(xd) + sqrt(N0*ofdm.fs/2)*randn(size(xd));

%% OFDM detection
% AdEq.W0 = exp(-1j*angle(F.H(ofdm.fc)));
[Xn, AGCn, W] = ofdm.detect(yk, AdEq, true);

%% Calculate BER
[ber_count, ~] = ofdm.count_ber([AdEq.Ntrain+sim.Ndiscard sim.Ndiscard]);

%% Gaussian approximation
Pnrx =  mean(abs(Xn).^2, 2)./abs(AGCn.*W).^2;
[ber_gauss, SNRndB] =  ofdm.estimate_ber(Pnrx.', Hch, varNoise, true);

ber_count 
ber_gauss

% profile off
% profile viewer
toc