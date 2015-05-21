%% Run transmission of dc-ofdm and check if with the selected clipping ratio 
%% is it possible to achieve target BER. Normally, when the error is higher
%% than the expected it's because the clippinig ratio is too small
% clear, clc, close all

format compact

addpath f           % functions path

rng('default')      % initiate default random number generator
rng('shuffle');     % Reinitialize the random number generator used by rand, randi, and randn with a seed based on the current time

%% Parameters to change
Nrept = 5;

sim.Nsymb = 2^14;                   % number of OFDM symbols

switch simulation_case
    case 1
        %% CASE 1: ENOB = 5, Preemphasis, CS = 16
        sim.type = 'preemphasis';                % Type of compensation at the transmitter {'preemphasis', 'palloc'}
        ofdm.CS = 16;                       % Constellation size (effective CS in case of variable bit loading)  
        sim.ENOB = 5;

        Fnl = (15:5:50)*1e9;

        RCLIPTX = [4.0 3.9 3.7 3.7 3.7 3.7 3.7 3.7];       % optmized clipping ratio
        RCLIPRX = [3.6 3.6 3.5 3.4 3.4 3.4 3.4 3.4];       % optmized clipping ratio
        
    case 2
        %% CASE 2: ENOB = 5, Palloc, CS = 16
        sim.type = 'palloc';                % Type of compensation at the transmitter {'preemphasis', 'palloc'}
        ofdm.CS = 16;                       % Constellation size (effective CS in case of variable bit loading)  
        sim.ENOB = 5;

        Fnl = (15:5:50)*1e9;

        RCLIPTX = [3.8 3.7 3.7 3.7 3.7 3.7 3.7 3.7];       % optmized clipping ratio
        RCLIPRX = [3.5 3.5 3.5 3.4 3.4 3.4 3.4 3.4];       % optmized clipping ratio
       
    case 3
        %% CASE 3: ENOB = 6, Preemphasis, CS = 16
        sim.type = 'preemphasis';                % Type of compensation at the transmitter {'preemphasis', 'palloc'}
        ofdm.CS = 16;                       % Constellation size (effective CS in case of variable bit loading)  
        sim.ENOB = 6;

        Fnl = (10:5:50)*1e9;

        RCLIPTX = [4.2 4.0 3.9 3.7 3.7 3.7 3.7 3.7 3.7];       % optmized clipping ratio
        RCLIPRX = [3.5 3.5 3.5 3.5 3.5 3.4 3.4 3.4 3.4];       % optmized clipping ratio
        
    case 4
        %% CASE 4: ENOB = 6, Preemphasis, CS = 64
        sim.type = 'preemphasis';                % Type of compensation at the transmitter {'preemphasis', 'palloc'}
        ofdm.CS = 64;                       % Constellation size (effective CS in case of variable bit loading)  
        sim.ENOB = 6;

        Fnl = (15:5:50)*1e9;

        RCLIPTX = [4.2 4.1 4.0 4.0 4.0 4.0 4.0 4.0];       % optmized clipping ratio
        RCLIPRX = [3.8 3.8 3.8 3.8 3.8 3.8 3.8 3.8];       % optmized clipping ratio
        
    case 5
        %% CASE 5: ENOB = 6, Palloc, CS = 16
        sim.type = 'palloc';                % Type of compensation at the transmitter {'preemphasis', 'palloc'}
        ofdm.CS = 16;                       % Constellation size (effective CS in case of variable bit loading)  
        sim.ENOB = 6;

        Fnl = (10:5:50)*1e9;

        RCLIPTX = [3.9 3.8 3.7 3.7 3.7 3.7 3.7 3.7 3.7];       % optmized clipping ratio
        RCLIPRX = [3.5 3.5 3.4 3.4 3.4 3.4 3.4 3.4 3.4];       % optmized clipping ratio

    case 6
        %% CASE 6:  ENOB = 6, Palloc, CS = 64
        sim.type = 'palloc';                % Type of compensation at the transmitter {'preemphasis', 'palloc'}
        ofdm.CS = 64;                       % Constellation size (effective CS in case of variable bit loading)  
        sim.ENOB = 6;

        Fnl = (10:5:50)*1e9;

        RCLIPTX = [4.4 4.3 4.1 4.1 3.9 3.9 3.9 3.9 3.9];       % optmized clipping ratio
        RCLIPRX = [3.5 3.5 3.4 3.4 3.4 3.4 3.4 3.4 3.4];       % optmized clipping ratio
        
    otherwise
        error('invalid option!')
end

%% Simulation
% Additional parameters
imp_resolution = 1e-9;  % Minimum resolution of the impulse response of the filters
                        % the last value of the truncated, normalized
                        % impulse response must be smaller than this value.
                        

verbose = false;

%%
%% Simulation parameters
sim.save = true;                   % save important data from simulation on file
sim.Pb = 1.8e-4;                    % Target BER
% sim.Nsymb = 2^14;                   % number of OFDM symbols
sim.Mct = 8;                       % oversampling rate for emulation of continuous time 
                                    % Note: Mct has to be at least 4 so DT filter design approximates CT well enough./
                                                                      
sim.quantiz = true;                % include quantization at both transmitter and receiver

%% Modulation parameters 
ofdm.ofdm = 'dc_ofdm';             % ofdm type
ofdm.full_dc = false;
ofdm.Nc = 64;                       % total number of subcarriers = FFT size (must be even)
ofdm.Nu = 52;                     % number of nonzero subcarriers (including complex conjugate)
ofdm.Rb = 107e9;                    % bit rate (total over two polarizations) (b/s)

ofdm.Ms = ofdm.Nc/ofdm.Nu;       % oversampling ratio
ofdm.Rs = 2*ofdm.Rb/log2(ofdm.CS);   % Symbol rate  
ofdm.B = ofdm.Nu/2*log2(ofdm.CS);    % Total number of bits per symbol

%% Transmitter parameters 
tx.kappa = 1;                         % current to optical power conversion (dc slope)

% Interpolation type. Could be one of the following
% {'ideal', 'butter', 'cheby1', 'ellipt', 'gaussian', 'bessel', 'fir'}
% Ideal is the ideal FIR interpolator design by interp
% The others refer to the type of filter that is used after ZOH
tx.filter = 'bessel';
tx.filter_order = 5;                % filter order
tx.filter_cutoff = 1/ofdm.Ms;       % Cut-off frequency normalized by ofdm.fs
tx.imp_length = 512;                % length of the impulse response of the filter

[bfilt, afilt] = design_filter(tx.filter, tx.filter_order, tx.filter_cutoff, sim.Mct);
                 
switch tx.filter
    case 'ideal'
        tx.gdac = bfilt/sum(bfilt);       
        tx.gdac_delay = grpdelay(bfilt, afilt, 1);
    otherwise
        bzoh = ones(1, sim.Mct);    % Grp delay = (Mct-1)/2
        bdac = conv(bfilt, bzoh);   % numerator of DAC transfer function (Gzoh x Gfilt) 
        adac = afilt;               % denominator of DAC transfer function (Gzoh x Gfilt) ZOH is FIR
        tx.gdac = impz(bdac, adac, tx.imp_length).';    

        tx.gdac = tx.gdac/sum(tx.gdac);
        tx.gdac_delay = grpdelay(bdac, adac, 1);    % calculates the group delay by approximating the filter as an FIR filter
                                                    % whose impulse response is given by tx.gdac   

        % Check if number of points used is high enough to attain desired
        % resolution
        assert(abs(tx.gdac(end)) < imp_resolution, 'tx.gdac length is not high enough');
end
                                             
%% Receiver parameters
rx.R = 1;                           % responsivity
rx.NEP = 30e-12;                    % Noise equivalent power of the TIA at the receiver (A/sqrt(Hz))
rx.Sth = rx.R^2*rx.NEP^2/2;            % two-sided psd of thermal noise at the receiver (Sth = N0/2)

% Antialiasing filter
rx.filter = 'gaussian';
rx.filter_order = 4;                            % filter order
rx.filter_cutoff = 1/ofdm.Ms;                   % Cut-off frequency normalized by ofdm.fs
rx.imp_length = 300;                            % length of the impulse response of the filter

[bfilt, afilt] = design_filter(rx.filter, rx.filter_order, rx.filter_cutoff, sim.Mct);

rx.gadc = impz(bfilt, afilt, rx.imp_length).';  % impulse response
rx.gadc = rx.gadc/sum(rx.gadc);                 % normalize so frequency response at 0 Hz is 0 dB
rx.gadc_delay = grpdelay(bfilt, afilt, 1);

% Check if number of points used is high enough to attain desired
% resolution
assert(abs(rx.gadc(end)) < imp_resolution, 'rx.gadc length is not high enough');

%% Iterate simulation varying the cut-off frequency of the 2nd-order filter of the transmitter
% Note: cut-off frequency changes the cyclic prefix length and thus chnages
% the sampling rate. 

% Initiliaze variables
PtxdBm = zeros(size(Fnl));
PtxedBm = zeros(size(Fnl));
bercount = zeros(size(Fnl));
berest = zeros(size(Fnl));
berint = zeros(length(Fnl), 2);

PP = zeros(length(Fnl), 2);

fprintf('------ %s, CS = %d, ENOB = %d ------\n', sim.type, ofdm.CS, sim.ENOB)

for k = 1:length(Fnl)
    fprintf('-------------- fnl = %d GHz --------------\n', Fnl(k)/1e9)
    sim.rcliptx = RCLIPTX(k);
    sim.rcliprx = RCLIPRX(k);
    
    % Transmitter parameters
    tx.fnl = Fnl(k);
    tx.Hl = @(f) 1./(1 + 2*1j*f/tx.fnl - (f/tx.fnl).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
    tx.hl = @(t) (2*pi*tx.fnl)^2*t(t >= 0).*exp(-2*pi*tx.fnl*t(t >= 0));
    tx.hl_delay = 2/(2*pi*tx.fnl);          % group delay of second-order filter in seconds

    % Calculate cyclic prefix
    [ofdm.Npre_os, ofdm.Nneg_os, ofdm.Npos_os] = cyclic_prefix(ofdm, tx, rx, sim);

    % Time and frequency scales
    ofdm.fs = ofdm.Rs*(ofdm.Nc + ofdm.Npre_os)/ofdm.Nu;               % sampling rate (Hz)
    ofdm.fsct = sim.Mct*ofdm.fs;                                      % sampling frequency to emulate continuous time (Hz)
    ofdm.fc = ofdm.fs/ofdm.Nc*(1:ofdm.Nu/2);                          % frequency at which subcarriers are located       
    
    dt = 1/ofdm.fsct;                                                 % time increment in emulating continuous time (s)
    Ntot = sim.Mct*(ofdm.Nc + ofdm.Npre_os)*sim.Nsymb;                % total number of points simulated in continuous time
    df = ofdm.fsct/Ntot;                                              % frequency increment in continuous time (Hz)                    

    sim.t = dt*(0:Ntot-1);                                            % continuous time scale
    sim.f = -ofdm.fsct/2:df:ofdm.fsct/2-df;                           % frequency range in continuous time                    
                               
      
    %%
    berc = zeros(1, Nrept);
    bere = zeros(1, Nrept);
    Ptxm = zeros(1, Nrept);
    Ptxe = zeros(1, Nrept);
    for kk = 1:Nrept
        [ber, P] = dc_ofdm(ofdm, tx, rx, sim, verbose);
   
        berc(kk) = ber.count;
        bere(kk) = ber.est;
        Ptxm(kk) = P.Ptx; 
        Ptxe(kk) = P.Ptxest; 
    end
    
    % format some important variables
    PtxdBm(k) = 10*log10(mean(Ptxm)/1e-3);
    PtxedBm(k) = 10*log10(mean(Ptxe)/1e-3);
    bercount(k) = log10(mean(berc));
    berest(k) = log10(mean(bere));
end

%% Required power for ACO-OFDM AWGN (with ideal filters, without CP, without dc bias)
snr = fzero(@(x) berqam(ofdm.CS, x) - sim.Pb, 20);      % snr to achieve target BER;
snrl = 10^(snr/10);

K = 1 - 2*qfunc(mean(RCLIPTX));

Pnrx = snrl*ofdm.Ms*ofdm.Rs*rx.Sth/ofdm.Nc;
Pnawgn = K.^2*Pnrx./abs(rx.R*tx.kappa)^2*ones(1, ofdm.Nu/2);
Pawgn = tx.kappa*sqrt(2*sum(Pnawgn))*mean(RCLIPTX);
PawgndBm = 10*log10(Pawgn/1e-3);

%% Required power for OOK AWGN
PookdBm = 10*log10(1/rx.R*qfuncinv(sim.Pb)*sqrt(rx.Sth*ofdm.Rb)/1e-3);

% Same as solving for square constellations
% qfuncinv((sim.Pb*log2(ofdm.CS)/4)/(1-1/sqrt(ofdm.CS)))/sqrt(3*pi*log2(ofdm.CS)/((ofdm.CS-1)*4*rx.Sth*ofdm.Rb))

% Power penalty with respect to NRZ-OOK assuming no bandwidth limitation.
% Oversampling penalty is included
PP_ofdm_dc = @ (M, r) 10*log10(r*sqrt(2*(M-1)./(3*log2(M)))); 

power_pen_ook_m = PtxdBm - PookdBm;
power_pen_ook_e = PtxedBm - PookdBm;

%% Figures
figure
subplot(211), hold on
plot(Fnl/1e9, bercount, '-sk', Fnl/1e9, berest, '-or')
plot(Fnl/1e9, log10(sim.Pb)*ones(size(Fnl)), 'k')
xlabel('Cut-off frequency (GHz)')
ylabel('log_{10}(BER)')
legend('BER counted', 'BER estimaded', 'Target BER')
axis([Fnl(1)/1e9 Fnl(end)/1e9 -4 -3])

subplot(212), hold on
plot(Fnl/1e9, power_pen_ook_m, '-sk', 'LineWidth', 2)
plot(Fnl/1e9, power_pen_ook_e, '-xb', 'LineWidth', 2)
plot(Fnl/1e9, (PawgndBm - PookdBm)*ones(size(Fnl)), '-r', 'LineWidth', 2)
xlabel('Cut-off frequency (GHz)', 'FontSize', 12)
ylabel('Optical Power Penalty (dB)', 'FontSize', 12)
legend('Measured', 'Estimated', 'AWGN')

if sim.save
    filepath = ['results/' ofdm.ofdm '/'];
    filename = [sim.type '_CS=' num2str(ofdm.CS) '_ENOB=' num2str(sim.ENOB)];
    
    saveas(gca, [filepath filename], 'fig');
end