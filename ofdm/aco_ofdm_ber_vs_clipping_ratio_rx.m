%% BER as a function of the clipping ratio at the receiver
clear, clc, close all

format compact

addpath f/           % functions path
rng('default')      % initiate default random number generator

%% Parameters to change
sim.ENOB = 5;   

cases = 3;          % select simulation

sim.quantiz = true;                % include quantization at both transmitter and receiver

sim.Nsymb = 2^14;

switch cases
    case 1 % Preemphasis, CS = 16
        sim.type = 'preemphasis';
        
        ofdm.CS = 16;
        
        Fnl = (10:5:50)*1e9;            % Modulator cutoff frequency
        
        RCLIPTX = [3.9 3.9 3.8 3.8 3.7 3.7 3.7 3.7 3.7]; % Optimized clipping ratio at the transmitter for each modulator frequency
        
        RCLIPRX = 3:0.5:6; % range to swipe
        
    case 2 % Preemphasis, CS = 64
        sim.type = 'preemphasis';
        
        ofdm.CS = 64;
        
        Fnl = (10:5:50)*1e9;
        
        RCLIPTX = [4.1 4.0 4.0 4.0 4.0 4.0 4.0 4.0 4.0];
        
        RCLIPRX = 3:0.2:4; % range to swipe
        
    case 3 % Palloc, CS = 16
        sim.type = 'palloc';
        
        ofdm.CS = 16;
        
        Fnl = (10:5:50)*1e9;
        
        RCLIPTX = [3.9 3.9 3.8 3.8 3.8 3.8 3.8 3.8 3.8];
        
        RCLIPRX = 2:0.25:3.5; %3:0.2:4; % range to swipe       
        
    case 4 % Palloc, CS = 64
        sim.type = 'palloc';
        
        ofdm.CS = 64;
        
        Fnl = 30e9; %(10:5:50)*1e9;
  
        RCLIPTX = [4.5 4.3 4.3 4.3 4.3 4.2 4.2 4.2 4.2];
        
        RCLIPRX = 3:0.2:4; % range to swipe 
        
    otherwise
        error('invalid option');
end

verbose = false;

Nave = 1; %Number of noise realizations

tol = 1e-2; % tolerance for the derivative of the BER. If the derivative becomes smaller than this value,
            % then the BER is assumed to be invariant with the clipping
            % ratio.

% Additional parameters
imp_resolution = 1e-9;  % Minimum resolution of the impulse response of the filters
                        % the last value of the truncated, normalized
                        % impulse response must be smaller than this value.

%%
%% Simulation parameters
sim.save = false;                   % save important data from simulation on file
% sim.type = 'preemphasis';                % Type of compensation at the transmitter {'preemphasis', 'palloc'}
sim.Pb = 1.8e-4;                    % Target BER
% sim.Nsymb = 2^14;                   % number of OFDM symbols
sim.Mct = 8;                        % oversampling rate for emulation of continuous time 
                                    % Note: Mct has to be at least 4 so DT filter design approximates CT well enough.

% Note: sim.rcliprx is only used if quantiz is set to true. Otherwise, it
% assumes ideal receiver (without quantization and clipping)
                                    
% sim.quantiz = false;                % include quantization at both transmitter and receiver
% sim.ENOB = 6;                       % effective number of bits
% sim.rcliprx = 6;                    % Clipping ratio of ACO-OFDM


%% Modulation parameters 
ofdm.ofdm = 'aco_ofdm';             % ofdm type
ofdm.Nc = 64;                       % total number of subcarriers = FFT size (must be even)
ofdm.Nu = 52/2;                     % number of nonzero subcarriers (including complex conjugate)
% ofdm.CS = 64;                       % Constellation size (effective CS in case of variable bit loading)    
ofdm.Rb = 107e9;                    % bit rate (total over two polarizations) (b/s)

ofdm.Ms = ofdm.Nc/(2*ofdm.Nu);       % oversampling ratio
ofdm.Rs = 4*ofdm.Rb/log2(ofdm.CS);   % Symbol rate  
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
rx.imp_length = 512;                            % length of the impulse response of the filter

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
PtxdBm = zeros(length(Fnl), length(RCLIPRX));
bercount = zeros(length(Fnl), length(RCLIPRX));
berest = zeros(length(Fnl), length(RCLIPRX));

for nn = 2:3 %length(Fnl)
    fprintf('-------------- fnl = %d GHz --------------\n', Fnl(nn)/1e9)
    
    tx.fnl = Fnl(nn);
    sim.rcliptx = RCLIPTX(nn);
    
    % Transmitter parameters
    tx.Hl = @(f) 1./(1 + 2*1j*f/tx.fnl - (f/tx.fnl).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
    tx.hl = @(t) (2*pi*tx.fnl)^2*t(t >= 0).*exp(-2*pi*tx.fnl*t(t >= 0));
    tx.hl_delay = 2/(2*pi*tx.fnl);          % group delay of second-order filter in seconds
    
    % Generate time and frequency scales
    % Calculate cyclic prefix
    [ofdm.Npre_os, ofdm.Nneg_os, ofdm.Npos_os] = cyclic_prefix(ofdm, tx, rx, sim);

    % Time and frequency scales
    ofdm.fs = ofdm.Rs*(ofdm.Nc + ofdm.Npre_os)/(2*ofdm.Nu);           % sampling rate (Hz)
    ofdm.fsct = sim.Mct*ofdm.fs;                                      % sampling frequency to emulate continuous time (Hz)
    ofdm.fc = ofdm.fs/ofdm.Nc*(1:2:ofdm.Nu);                          % frequency at which subcarriers are located       

    dt = 1/ofdm.fsct;                                                 % time increment in emulating continuous time (s)
    Ntot = sim.Mct*(ofdm.Nc + ofdm.Npre_os)*sim.Nsymb;                % total number of points simulated in continuous time
    df = ofdm.fsct/Ntot;                                              % frequency increment in continuous time (Hz)                    

    sim.t = dt*(0:Ntot-1);                                            % continuous time scale
    sim.f = -ofdm.fsct/2:df:ofdm.fsct/2-df;                           % frequency range in continuous time      

    % reset and save rng seed
    rng('shuffle');
    rng_seed = rng;
    
    for k = 1:length(RCLIPRX)
        fprintf('--- r = %.2f --- \n', RCLIPRX(k))
        sim.rcliprx = RCLIPRX(k);                    % Clipping ratio of ACO-OFDM

        PtxdBmv = zeros(1, Nave);
        berc = zeros(1, Nave);
        bere = zeros(1, Nave);

        %% 
        rng(rng_seed); % make sure that the noise and ofdm symbol is consistent for all clipping ratios tested
        for kk = 1:Nave
            try
                [ber, P] = aco_ofdm(ofdm, tx, rx, sim, verbose);
            catch err
                disp(err.message)
                fprintf('fnl = %d GHz, r = %.2f\n', tx.fnl/1e9, sim.rcliprx);
                break;
            end
            
            PtxdBmv(kk) = 10*log10(P.Ptx/1e-3);
            berc(kk) = log10(ber.count);
            bere(kk) = log10(ber.est);
        end

        % format some important variables
        PtxdBm(nn, k) = mean(PtxdBmv);
        bercount(nn, k) = mean(berc);
        berest(nn, k) = mean(bere);
        % each row corresponds to one cutoff frequency and each column
        % corresponds to one clipping ratio
    end

    %% Required power for ACO-OFDM AWGN (with ideal filters, without CP, without dc bias)
    snr = fzero(@(x) berqam(ofdm.CS, x) - sim.Pb, 20);      % snr to achieve target BER;
    snrl = 10^(snr/10);

    Pnrx = snrl*ofdm.Ms*ofdm.Rs*rx.Sth/ofdm.Nc;
    Pnawgn = 4*Pnrx./abs(rx.R*tx.kappa)^2*ones(1, ofdm.Nu/2);
    Pawgn = tx.kappa*1/sqrt(2*pi)*sqrt(2*sum(Pnawgn)); 
    PawgndBm = 10*log10(Pawgn/1e-3);

    %% Required power for OOK AWGN
    PookdBm = 10*log10(1/rx.R*qfuncinv(sim.Pb)*sqrt(rx.Sth*ofdm.Rb)/1e-3);

    % Same as solving for square constellations
    % qfuncinv((sim.Pb*log2(ofdm.CS)/4)/(1-1/sqrt(ofdm.CS)))/sqrt(3*pi*log2(ofdm.CS)/((ofdm.CS-1)*4*rx.Sth*ofdm.Rb))

    % Power penalty with respect to NRZ-OOK assuming no bandwidth limitation.
    % Oversampling penalty is included
    PP_ofdm_aco = @ (M) 10*log10(sqrt(4*(M-1)./(3*pi*log2(M)))); 

    power_pen_ook = PtxdBm - PookdBm;

    %% Estimate optimal clipping ratio 
    % if no clipping ratio value makes it work for this frequency, 
    % then move to the next frequency
    if all(bercount(nn,:) == 0)
        fprintf('----- There is no clipping ratio that makes it work at fnl = %d GHz -----\n', tx.fnl/1e9);
        continue
    elseif length(bercount(nn, bercount(nn,:) ~= 0)) >= 2
        nz_pos = (bercount(nn,:) ~= 0);
        RCLIPRX_nz = RCLIPRX(bercount(nn,:) ~= 0);
        rrx = linspace(RCLIPRX_nz(1), RCLIPRX_nz(end));
        bercount_int = spline(RCLIPRX_nz, bercount(nn, bercount(nn,:) ~= 0), rrx);
        drrx = abs(rrx(1)-rrx(2));
        dber = diff(bercount_int)/drrx;
        
        %% find optimal when derivative of the curve goes below a certain threshold for the first time
        if isempty(find(abs(dber) < tol, 1, 'first'))
            rrxopt(nn) = NaN; % put not-a-number if doesn't exist optimal
        else
            rrxopt(nn) = rrx(find(abs(dber) < tol, 1, 'first'));
        end
        
        %% find optimal assuming that derivative of the curve stays below a certain threshold
        %     loc = (abs(dber) < tol);
        % 
        %     % optimal clipping ratio is the one for each all the following dber are
        %     % less the tolerance value.
        %     for k = 1:length(loc)
        %         if all(loc(k:end))
        %             rrxopt(nn) = rrx(k);
        %             break;
        %         else
        %             continue;
        %         end
        %     end
    else
        rrxopt(nn) = bercount(nn, bercount(nn,:) ~= 0);
    end
    
    %% Figures
    figure
    subplot(211), hold on
    plot(RCLIPRX_nz, bercount(nn, nz_pos), 'sk', RCLIPRX_nz, berest(nn,nz_pos), '-or')
    plot(RCLIPRX, log10(sim.Pb)*ones(size(RCLIPRX)), 'k')
    plot(rrx, bercount_int, '-k')
    if ~isnan(rrxopt(nn)) && rrxopt(nn) ~= 0
        plot(rrxopt(nn), bercount_int(rrx == rrxopt(nn)), '*k', 'MarkerSize', 10);
    end

    xlabel('Clipping ratio', 'FontSize', 12)
    ylabel('log_{10}(BER)', 'FontSize', 12)
    legend('BER counted', 'BER estimaded', 'Target BER')
    axis([RCLIPRX([1 end]) -4 -3])
    grid on
    title(sprintf('fnl = %d GHz', tx.fnl/1e9));

    subplot(212), hold on
    plot(rrx(1:end-1), abs(dber), 'k')
    plot(rrx([1 end-1]), [tol tol], '-r')
    xlabel('Clipping ratio', 'FontSize', 12)
    ylabel('|dlog_{10}(BER)/dr|', 'FontSize', 12)
    axis([RCLIPRX([1 end]) 0 2*tol])
    
    
    if sim.save
        filepath = ['results/' ofdm.ofdm '/'];
        if sim.quantiz
            filename = ['clip_opt_' ofdm.ofdm '_' sim.type '_CS=' num2str(ofdm.CS) '_ENOB=' num2str(sim.ENOB) '_Fnl=' num2str(Fnl(nn)/1e9) 'GHz'];
        else
            filename = ['clip_opt_' ofdm.ofdm '_' sim.type '_CS=' num2str(ofdm.CS) '_Fnl=' num2str(Fnl/1e9) 'GHz'];
        end
        
        % save plot
        saveas(gca, [filepath filename], 'png')
    end
    
end

if sim.save
    filepath = ['results/' ofdm.ofdm '/'];
    if sim.quantiz
        filename = ['clip_opt_' ofdm.ofdm '_' sim.type '_CS=' num2str(ofdm.CS) '_ENOB=' num2str(sim.ENOB)];
    else
        filename = ['clip_opt_' ofdm.ofdm '_' sim.type '_CS=' num2str(ofdm.CS)];
    end
    
    % remove heavy data before saving
    sim = rmfield(sim, {'f', 't'});

    % save data
    save([filepath filename])
end
