%% Generate and detect DC-OFDM (fiber effects are not included)
% input:
% ofdm (OFDM parameters)
% tx (transmitter parameters)
% fiber (fiber parameters)
% rx (receiver parameters)
% sim (simulation parameters)
% verbose = display figures?

% output:
% ber: counted (including 95%-confidence intervals), and estimated
% P: structure containing
% % Ptx = average optical fiber (measured)
% % Ptx_est = average optical fiber (estimated)
% % Pn = power allocation (power at each subcarrier)
% pp: struct of estimated power penalties: pp.{dcbias, cp, nonflat}
% SNRn: SNR at nth subcarrier
% CS: vector of constellation size for each subcarrier
% approx: struct containing checks for approximations made. 
% % Default usage: approx.variable_name = [theory, measured]

% Function assumes variable number of input arguments just to be compatible
% with previous codes
function [ber, P, pp, SNRn, CS, approx] = dc_ofdm(varargin)

if nargin == 5 % if 5 input elements, then fiber propagation is not simulated
    [ofdm, tx, rx, sim, verbose] = varargin{:};
    
    Hcd = 1;

elseif nargin == 6 % if 6 input arguments, then fiber struct was passed
    [ofdm, tx, fiber, rx, sim, verbose] = varargin{:};
    
    % Dispersion parameter to beta2 (lamb0 is the reference wavelength)
    c = 299792458;  % speed of light
    D2beta2 = @(D, lamb) -D*lamb^2/(2*pi*c); 
    
    fiber.beta2 = D2beta2(fiber.D(tx.lamb), tx.lamb);
    
    % CD frequency response
    theta = -1/2*fiber.beta2*(2*pi*ofdm.fc).^2*fiber.L; % theta = -1/2*beta2*w.^2*L
    
    Hcd = cos(theta) - tx.alpha*sin(theta);        % fiber small-signal frequency response
   
else
    error('Invalid number of input arguments.');
end

assert(strcmp(ofdm.ofdm, 'dc_ofdm'), 'Simulation was not set for dc-ofdm');

%% Define and redefine some variables
Ndis = 128;                       % Number of symbols to be discarded from results (both at the beginning and at the end)

f = sim.f;
Nc = ofdm.Nc; 
Nu = ofdm.Nu;
Npre_os = ofdm.Npre_os;
Npos_os = ofdm.Npos_os;
Nneg_os = ofdm.Nneg_os;
Nsymb = sim.Nsymb;

K = 1 - 2*qfunc(sim.rcliptx);                % Amplitude attenuation due to clipping = (1-Q(r))

% Remove group delay of modulator frequency response
Hl = tx.kappa*tx.Hl(f);
Hl = Hl.*exp(1j*2*pi*f*tx.hl_delay);

% The group delay is removed from the frequency response of the ADC and DAC
Gdac = fft(tx.gdac, length(f));
Gdac = fftshift(Gdac).*exp(1j*2*pi*f*tx.gdac_delay/ofdm.fsct);                        

Gadc = fft(rx.gadc, length(f));
Gadc = fftshift(Gadc).*exp(1j*2*pi*f*rx.gadc_delay/ofdm.fsct);   

% At the subcarriers frequencies
nfc = length(f)/2+1 + (1:Nu/2)*(Nc+Npre_os)/Nc*Nsymb;           % index for subcarrier frequencies  

tx.Gdac = Gdac(nfc); % DAC frequency response at the subcarriers frequencies
tx.Hmod = Hl(nfc); % Modulator frequency response at the subcarriers frequencies
fiber.Hcd = Hcd; % Fiber frequency response at the subcarriers frequencies
rx.Gadc = Gadc(nfc); % ADC frequency response at the subcarriers frequencies
Gch = K*tx.Gdac.*tx.kappa.*tx.Hmod.*fiber.Hcd.*rx.R.*rx.Gadc;            % Frequency response of the channel at the subcarriers

% Include attenuation
Gch = 10^(-fiber.att*fiber.L/1e4)*Gch;

%% Define how to calculate dc bias
% (i) full dc bias: sigtx*rcliptx
% (ii) reduced dc bias: sigtxdac*rcliptx
if isfield(ofdm, 'full_dc') == 0
    % add full dc-bias if ofdm.full_dc is true or if it was not defined (compatibility with old codes)
    warning('undefined variable ofdm.full_dc. Assuming full dc-bias simulation!');
        
    ofdm.full_dc = true;
end

if ofdm.full_dc
    ofdm.dc_bias = @(Pn) sim.rcliptx*sqrt(sum(2*Pn));
else
    ofdm.dc_bias = @(Pn) sim.rcliptx*sqrt(sum(2*Pn.*abs(Gdac(nfc)).^2));
end

%% Preemphasis at the transmitter
switch sim.type
    case 'preemphasis'
        [Pn, CS, aux] = preemphasis(ofdm, tx, fiber, rx, Gch, sim);
    %% Optimal power allocation and bit loading   
    case 'palloc'
        [Pn, CS, aux] = palloc(ofdm, tx, fiber, rx, Gch, sim);
    otherwise 
        disp('invalid option') 
end

% Signal std at transmitter and receiver
sigtx = sqrt(2*sum(Pn));
sigrx = sqrt(2*sum(Pn.*K^2.*abs(Gch).^2));

% Check if power is at reasonable levels
assert(sigtx < 1e-3, 'Power is too high: %.2G (mW)\nCannot improve performance by increasing power.\n', 1e3*sqrt(2*sum(Pn)))

Ptx_est = tx.kappa*ofdm.dc_bias(Pn);
Prx_est = Ptx_est/10^(fiber.att*fiber.L/1e4);

% Calculate the average power of a CS-QAM constellation with dmin = 2 (used
% to normalize symbol power after qammod)
Pqam = zeros(1, Nu/2);
for k = 1:Nu/2
    Pqam(k) = mean(abs(qammod(0:CS(k)-1, CS(k), 0, 'gray')).^2);
end

if ~isfield(sim, 'Navg')
    sim.Navg = 1; % if not defined averages only once
end

for navg = 1:sim.Navg

    %% Generate OFDM signal at chip rate (done in DSP)
    dataTX = zeros(Nu/2, Nsymb);
    dataTXm = zeros(Nu/2, Nsymb);
    for kk = 1:Nu/2
        dataTX(kk,:) = randi([0 CS(kk)-1], [1 Nsymb]);              % data to be encoded (symbols columnwise)
        dataTXm(kk,:) = qammod(dataTX(kk,:), CS(kk), 0, 'gray');    % encoded QAM symbols to be modulated onto subcarriers
        dataTXm(kk,:) = sqrt(Pn(kk))*dataTXm(kk,:)/sqrt(Pqam(kk));  % scale constellation so that Pn is the power at nth subcarrier (E(|Xn|^2) = Pn)
    end

    % zero-pad and ensure Hermitian symmetry
    % -f(N) ... -f(1) f(0) f(1) ... f(N-1) => (0 ... 0, X(-n)*, 0, X(n), 0 ... 0)
    Xn = ifftshift([zeros((Nc-Nu)/2, Nsymb); ...     % zero-pad neg. frequencies
                  flipud(conj(dataTXm));...            % data*(-n) (Nu) 
                  zeros(1, Nsymb); ...                % 0 at f == 0 (1)
                  dataTXm; ...                         % data(n)  (Nu)
                  zeros((Nc-Nu)/2 - 1, Nsymb)], 1);   % zero-pad pos. frequencies

    % Perform ifft (Nc-IDFT) columnwise
    xn = Nc*ifft(Xn, Nc, 1); 

    % Insert cyclic prefix
    xncp = [xn(end-Npre_os+1:end, :); xn]; % insert cyclic prefix

    % Parallel to serial
    xncp = reshape(xncp, 1, (Nc + Npre_os)*Nsymb); % time-domain ofdm signal w/ cyclic prefix

    %% Clipping
    % Note: both positive and negative tails are clipped
    xncpc = xncp;
    xncpc(xncp < -sim.rcliptx*sigtx) = -sim.rcliptx*sigtx;
    xncpc(xncp > sim.rcliptx*sigtx) = sim.rcliptx*sigtx;

    % Check approximations
    approx.clip_prob = 2*[qfunc(sim.rcliptx) mean(xncp > sim.rcliptx*sigtx)]; % Clipping probability


    %% Quantize at the transmitter
    if sim.quantiz
        % With quantization (clip -> quantiz -> interp)  
        % Quantiz
        yqtx = linspace(-sim.rcliptx*sigtx, sim.rcliptx*sigtx, 2^sim.ENOB);     % Quantization levels
        delta1 = abs(yqtx(2)-yqtx(1));              % sim.rcliptx*sigtx/(2^sim.ENOB-1);        % (2^sim.ENOB-1) because of clipping done before
        [~, xncpq, varQ1] = quantiz(xncpc, yqtx(1:end-1) + delta1/2, yqtx);

        QuantErrTX = xncpc - xncpq; % used to plot histogram

        % check approximations
        approx.delta1 = [aux.delta1th, delta1];
        approx.varQ1 = [aux.varQ1th, varQ1];
    else 
        % Without quantization (clip -> interp)
        % No quantization
        varQ1 = 0;
        xncpq = xncpc;
    end

    %% Interpolation
    % Sampling rate expansion
    xt = zeros(1, sim.Mct*sim.Nsymb*(Npre_os + Nc));
    xt(1:sim.Mct:end) = xncpq;

    % Filter (interpolator + ZOH)
    xt = real(ifft(ifftshift(sim.Mct*Gdac).*fft(xt)));
    % Note: multiply by sim.Mct because DAC filter should have gain sim.Mct
    % to compensate for 1/sim.Mct due sampling rate expansion

    %% Add dc bias (dc_bias determines whether it's added full dc bias or clipped dc bias)
    xt = xt + ofdm.dc_bias(Pn);

    % Ensure that signal is non-negative
    xtc = xt;
    xtc(xt < 0) = 0;

    %% Apply frequency response of the laser
    Pt = real(ifft(ifftshift(Hl).*fft(xtc)));  % Note: real() is used to remove residual imag
                                                % part that appears due to numerical error        

    %% Modulator Nonlinearity
    % To be implemented

    
    %% Add additional dc-bias due to finite exctinction ratio
    if isfield(tx, 'rexdB') % check if exctinction ratio (rex) was defined
        rex = 10^(tx.rexdB/10);
        
        Pmin_rex = 2*sim.rcliptx*sigtx/(rex - 1); %% additional dc bias
        
        Pt = Pt + Pmin_rex; 
                                                
    end
    
    % Calculate equivalent transmitted power
    Ptx  = mean(Pt);   % measured
    
    %% Adjust transmitted power
    if isfield(tx, 'Ptx') % if Ptx was provided, then scale signal to desired Ptx
        Pt = Pt*tx.Ptx/Ptx;
        Pn = (tx.Ptx/Ptx)^2*Pn;
        Ptx = tx.Ptx;
        warning('Average power readjusted')
    end

    %% Add intensity noise
    varrin = 0;
    if sim.RIN
        % Calculate RIN PSD, which depends on the instantaneous power
        Srin = 10^(tx.RIN/10)*Pt.^2;

        % noise. Srin is two-sided psd
        wrin = sqrt(Srin.*ofdm.fsct).*randn(size(Pt));
        
        % variance 
        varrin = var(wrin)/sim.Mct;

        Pt = Pt + wrin;
    end 

    %% Fiber propagation
    if exist('fiber', 'var') == 1 && fiber.D(tx.lamb)*fiber.L ~= 0
        disp('Including chromatic dispersion.')

        % Calculate electric field including chirp
        dphi = tx.alpha/2*log(Pt);         % Phase variation due to chirp (only transient chirp is considered)

        Et = sqrt(Pt).*exp(1j*dphi);

        % Electric field after fiber propagation
        Dw = -1j*D2beta2(fiber.D(tx.lamb), tx.lamb)/2*(2*pi*ifftshift(sim.f)).^2;
        Etz = ifft(exp(fiber.L*Dw).*fft(Et));

        % Received power 
        Ptz = 10^(-fiber.att*fiber.L/1e4)*abs(Etz).^2;
    else
        Ptz = Pt;
    end


    %% Direct detection
    Its = rx.R*Ptz;

    %% Add thermal noise
    wt = sqrt(rx.Sth*ofdm.fsct)*randn(1, sim.Mct*sim.Nsymb*(Npre_os + Nc));
    varthermal = var(wt)/sim.Mct;
    It = Its + wt;

    %% Add shot noise
    varshot = 0;
    if sim.shot
        q = 1.60217657e-19;      % electron charge (C)
        Id = 0;                  % dark current

        % Instataneous received power considering only attenuation from the fiber   
        Sshot = 2*q*(rx.R*Ptz + Id);     % one-sided shot noise PSD

        % Frequency is divided by two because PSD is one-sided
        ws = sqrt(Sshot*ofdm.fsct/2).*randn(size(It));
        
        % variance of shot noise at signal rate
        varshot = var(ws)/sim.Mct;
        
        wt = wt + ws;
        It = It + ws;
    end

    % Antialiasing filter of the ADC
    It = real(ifft(fft(It).*ifftshift(Gadc))); % signal + noise
    Its = real(ifft(fft(Its).*ifftshift(Gadc))); % signal only
    wt = real(ifft(fft(wt).*ifftshift(Gadc))); % noise only
    % Note: we will treat noise separately as it makes it easier to estimate
    % SNR (etc). If all the operations that follow are linear this should make
    % no difference

    %% Back to discrete time
    % Resample signal once per chip        
    yncp = It(1:sim.Mct:end);       % yn with cyclic prefix without noise  
    yncps = Its(1:sim.Mct:end);     % yn with cyclic prefix without noise     
    wncp = wt(1:sim.Mct:end);

    %% Quantize at the receiver
    if sim.quantiz    
        % remove dc bias
        yncp = yncp - mean(yncp);

        % quantization levels
        yqrx = linspace(-sim.rcliprx*sigrx, sim.rcliprx*sigrx, 2^sim.ENOB); % Quantization levels
        delta2 = abs(yqrx(2) - yqrx(1)); 

        % clip before quantization
        yncpc = yncp;
        yncpc(yncp < yqrx(1)) = yqrx(1);
        yncpc(yncp > yqrx(end)) = yqrx(end);

        % Quantize
        [~, yncpq, varQ2] = quantiz(yncpc, yqrx(1:end-1) + delta2/2, yqrx);

        QuantErrRX = yncpc - yncpq; % used to plot histogram if verbose = true

        % Check approximations
        approx.delta2 = [aux.delta2th, delta2];
        approx.varQ2 = [aux.varQ2th, varQ2];
    else
        % No quantization
        varQ2 = 0;
        yncpq = yncp - mean(yncp);
    end

    %%

    % reshape into matrix form
    yncpq = reshape(yncpq, Nc + Npre_os, Nsymb);      % signal + noise
    yncps = reshape(yncps, Nc + Npre_os, Nsymb);      % signal
    wncp = reshape(wncp, Nc + Npre_os, Nsymb);        % noise

    % Remove cyclic prefix
    yn = circshift(yncpq(Npos_os+1:end-Nneg_os, :), -Nneg_os);  % signal + noise
    yns = circshift(yncps(Npos_os+1:end-Nneg_os, :), -Nneg_os); % signal
    wn = circshift(wncp(Npos_os+1:end-Nneg_os, :), -Nneg_os);   % noise

    % Demodulate symbols 1/N*DFT()
    Yn = fft(yn, Nc, 1)/Nc;             % signal + noise      
    Yns = fft(yns, Nc, 1)/Nc;           % signal 
    Wn = fft(wn, Nc, 1)/Nc;             % noise

    % Calculate SNR at each subcarrier
    Yns = mean(abs(Yns).^2, 2);
    Yns = Yns(1 + (1:Nu/2)).';      % Only consider positive-frequency subcarriers

    Wn = mean(abs(Wn).^2, 2); 
    Wn = Wn(1 + (1:Nu/2)).';        % Only consider positive-frequency subcarriers

    SNRn = 10*log10(Yns./(1/Nc*abs(rx.Gadc).^2.*(varthermal + varshot + varrin*abs(10^(-fiber.att*fiber.L/1e4)*fiber.Hcd.*rx.R).^2) + varQ1/Nc*abs(Gch).^2 + varQ2/Nc)); % Measured SNR

    % used subcarrier amplitudes (complex conjugate subcarriers are ignored)
    dataRXm = Yn(1 + (1:Nu/2), Ndis+1:end-Ndis);            

    % AGC caculates what should be the scaling factor to normalize all
    % constellations so that dmin = 2 (value expected by qamdemod)
    AGCn = sqrt(Pqam)./(K*sqrt(Pn).*Gch);
    % Note: Factor of K appears due to clipping

    dataRX = zeros(Nu/2, Nsymb-2*Ndis);
    numerr = zeros(1, Nu/2);
    bn = log2(CS);
    for kk = 1:Nu/2
    %     scatterplot(dataRXm(kk,:))
        dataRXm(kk,:) = AGCn(kk)*dataRXm(kk,:);
    %     scatterplot(dataRXm(kk,:))
        dataRX(kk, :) = qamdemod(dataRXm(kk,:), CS(kk), 0, 'gray');         % encoded QAM symbols to be modulated onto subcarriers
        numerr(kk) = biterr(dataTX(kk, Ndis+1:end-Ndis), dataRX(kk,:), bn(kk));
    end

    % average BER
    ber.count(navg) = sum(numerr)/(sim.Nsymb*sum(bn));

    % 95% confidence intervals for the counted BER
    [~, interval] = berconfint(numerr, sim.Nsymb*bn);
    ber.interval(navg, :) = [min(interval(:,1)), max(interval(:,2))];

    % Estimate BER from the SNR at each subcarrier
    % Note that to estimate the BER we use SNRn which is calculated from the
    % data measured in the simulation. If we had used SNRnest (calculated from
    % the theoretical values) the results would match the target BER perfectly
    berest = zeros(size(SNRn));
    for kk = 1:length(SNRn)
        berest(kk) = berqam(CS(kk), SNRn(kk));
    end

    ber.est(navg) = sum(berest.*bn)/sum(bn);
end

ber.count = mean(ber.count);
ber.est = mean(ber.est);
ber.interval = mean(ber.interval, 2);

%% Power
P.Ptx = Ptx;
P.Pn = Pn;
P.Ptx_est = Ptx_est;

%% Approximations
approx.sigtx = [sigtx std(xncp)];
approx.sigrx = [sigrx std(yncp)];
approx.Ptx = [Ptx_est Ptx];

%% Penalties
% Penalty due to nonflat frequency response of the channel (DAC * Hl * ADC)
pp.nonflat = 10*log10(1/sqrt(2*pi)*sqrt(2*sum(K^2./abs(Gch).^2))); 
% Note: the factor K^2 appears multiplying rather than dividing the power
% because the average power depends is given by sigma/(2*pi), where sigma
% is the power before clipping.

% Cyclic prefix
pp.cp = 10*log10(sqrt((Nc+ofdm.Npre_os)/Nc));

%% Show plots if user set verbose to true
if verbose 
    disp('Approximations: [theory, measured]')
    approx
    
    plot_stuff
end

