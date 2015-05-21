%% Generate and detect ACO-OFDM (fiber effects are not included)
% input:
% ofdm (OFDM parameters)
% tx (transmitter)
% rx (receiver)
% sim (simulation parameters)
% verbose = display figures?

% output:
% ber: counted (including 95%-confidence intervals), and estimated
% Ptx: average optical fiber
% pp: struct of estimated power penalties: pp.{dcbias, cp, nonflat}
% SNRn: SNR at nth subcarrier
% Pn: power at nth subcarrier
% CS: vector of constellation size for each subcarrier
% approx: struct containing checks for approximations made. 
% Default usage: approx.variable_name = [theory, measured]
function [ber, P, pp, SNRn, CS, approx] = aco_ofdm(ofdm, tx, rx, sim, verbose)

assert(strcmp(ofdm.ofdm, 'aco_ofdm'), 'Simulation was not set for aco-ofdm');

%% Define and redefine some variables
Ndis = 128;                       % Number of symbols to be discarded from results (both at the beginning and at the end)
Nit_piterative = 50;            % Number of iterations used to estimate required power at subcarriers to achieve target BER for the case of quantization.

f = sim.f;
Nc = ofdm.Nc; 
Nu = ofdm.Nu;
Npre_os = ofdm.Npre_os;
Npos_os = ofdm.Npos_os;
Nneg_os = ofdm.Nneg_os;
Nsymb = sim.Nsymb;

K = 1/2 - qfunc(sim.rcliptx);                % Amplitude attenuation due to clipping = (1-Q(r=0))

% Remove group delay from modulator frequency response
Hl = tx.kappa*tx.Hl(f);
Hl = Hl.*exp(1j*2*pi*f*tx.hl_delay);

% Remove group delay from ADC and DAC
Gdac = fft(tx.gdac, length(f));
Gdac = fftshift(Gdac).*exp(1j*2*pi*f*tx.gdac_delay/ofdm.fsct);                        

Gadc = fft(rx.gadc, length(f));
Gadc = fftshift(Gadc).*exp(1j*2*pi*f*rx.gadc_delay/ofdm.fsct);   

% At the subcarriers frequencies
nfc = length(f)/2+1 + (1:2:Nu)*(Nc+Npre_os)/Nc*Nsymb;           % index for subcarrier frequencies          
Gch = Gdac(nfc).*tx.kappa.*Hl(nfc).*rx.R.*Gadc(nfc);            % Frequency response of the channel at the subcarriers

%% Preemphasis at the transmitter
switch sim.type
    case 'preemphasis'
        CS = ofdm.CS*ones(1,Nu/2);

        % SNR to achieve target BER sim.Pb
        snrdB = fzero(@(x) berqam(ofdm.CS, x) - sim.Pb, 20);
        snrl = 10^(snrdB/10);

        if sim.quantiz
            % Iterate to get correct power. For kk = 1, it assumes that the
            % power is zero i.e., there is no noise due to quantization
            varQ1th = 0;
            varQ2th = 0;
            for kk = 1:Nit_piterative
                Pnrx = snrl*(ofdm.fs*rx.Sth*abs(Gadc(nfc)).^2 + varQ1th*abs(Gch).^2 + varQ2th)/Nc;                                               

                Pn = Pnrx./(K^2*abs(Gch).^2);
                
                % Signal std at transmitter and receiver
                sigtx = sqrt(2*sum(Pn));
                sigrx = sqrt(2*sum(Pn.*abs(Gch).^2));
                % Note: in aco-ofdm the value sigrx does not correspond to
                % the std of the signal at the receiver. Instead, it's an
                % approximation to the right tail of the pdf of the signal.
 
                % Quantization noise variance
                delta1th = sim.rcliptx*sigtx/(2^sim.ENOB-1);
                delta2th = (sigtx/sqrt(2*pi) + sim.rcliprx*sigrx)/(2^sim.ENOB-1); 
                varQ1th = (1/2 - qfunc(sim.rcliptx))*delta1th^2/12;     % quantization at the transmitter
                varQ2th = (1 - qfunc(sim.rcliprx))*delta2th^2/12;       % quantization at the receiver
            end
        else
            % Calculate equivalent power at each subcarrier at the receiver to
            % achieve desired SNR
            Pnrx = snrl*ofdm.fs*rx.Sth*abs(Gadc(nfc)).^2/Nc;
            % Power at each subcarrier at the transmitter
            Pn = Pnrx./abs(K*Gch).^2;
            
            % No quantization
            varQ1th = 0;
            varQ2th = 0;
        end

    %% Optimal power allocation and bit loading   
    case 'palloc'
        if sim.quantiz
            varQ1th = 0;
            varQ2th = 0;
            for kk = 1:Nit_piterative
                % Gain-to-noise ratio
                GNR = K^2*Nc*abs(Gch).^2./(ofdm.fs*rx.Sth*abs(Gadc(nfc)).^2 + varQ1th*abs(Gch).^2 + varQ2th); 

                % Run Levin-Campello algorithm to determine optimal power allocation
                [~, CS, Pn] = Levin_Campello_MA(ofdm.B, 1, sim.Pb, GNR); 

                % Signal std at transmitter and receiver
                sigtx = sqrt(2*sum(Pn));
                sigrx = sqrt(2*sum(Pn.*abs(Gch).^2));
                % Note: in aco-ofdm the value sigrx does not correspond to
                % the std of the signal at the receiver. Instead, it's an
                % approximation to the right tail of the pdf of the signal.
                
                % Quantization noise variance
                delta1th = sim.rcliptx*sigtx/(2^sim.ENOB-1);
                delta2th = (sigtx/sqrt(2*pi) + sim.rcliprx*sigrx)/(2^sim.ENOB-1); 
                varQ1th = (1/2 - qfunc(sim.rcliptx))*delta1th^2/12;     % quantization at the transmitter
                varQ2th = (1 - qfunc(sim.rcliprx))*delta2th^2/12;       % quantization at the receiver
            end
        else
            % Gain-to-noise ratio
            GNR = K^2*Nc*abs(Gch).^2./(ofdm.fs*rx.Sth*abs(Gadc(nfc)).^2); 

            % Run Levin-Campello algorithm to determine optimal power allocation
            [~, CS, Pn] = Levin_Campello_MA(ofdm.B, 1, sim.Pb, GNR);

            % No quantization
            varQ1th = 0;
            varQ2th = 0;
        end 
    otherwise 
        disp('invalid option') 
end

% Signal std at transmitter and receiver
sigtx = sqrt(2*sum(Pn));
sigrx = sqrt(2*sum(Pn.*abs(Gch).^2));
% Note: in aco-ofdm the value sigrx does not correspond to
% the std of the signal at the receiver. Instead, it's an
% approximation to the right tail of the pdf of the signal.

% Estimated avarage power
Ptxest = 1/sqrt(2*pi)*sqrt(2*sum(Pn));

% Check if power is at reasonable levels
assert(Ptxest < 10e-3, 'Power is too high: %.2G (mW)\nCannot improve performance by increasing power.\n', 1e3*Ptxest)
    
% Calculate the average power of a CS-QAM constellation with dmin = 2 (used
% to normalize symbol power after qammod)
Pqam = zeros(1, Nu/2);
for k = 1:Nu/2
    Pqam(k) = mean(abs(qammod(0:CS(k)-1, CS(k), 0, 'gray')).^2);
end

%% Generate OFDM signal at chip rate (done in DSP)
dataTX = zeros(Nu/2, Nsymb);
dataTXm = zeros(Nu/2, Nsymb);
for kk = 1:Nu/2
    dataTX(kk,:) = randi([0 CS(kk)-1], [1 Nsymb]);              % data to be encoded (symbols columnwise)
    dataTXm(kk,:) = qammod(dataTX(kk,:), CS(kk), 0, 'gray');    % encoded QAM symbols to be modulated onto subcarriers
    dataTXm(kk,:) = sqrt(Pn(kk))*dataTXm(kk,:)/sqrt(Pqam(kk));  % scale constellation so that Pn is the power at nth subcarrier (E(|Xn|^2) = Pn)
end

% Set even subcarriers to zero and assign data to odd subcarriers
dataTXm2 = zeros(Nu, Nsymb);
dataTXm2(1:2:end, :) = dataTXm;

% zero-pad and ensure Hermitian symmetry
% -f(N) ... -f(1) f(0) f(1) ... f(N-1) => (0 ... 0, X(-n)*, 0, X(n), 0 ... 0)
Xn = ifftshift([zeros((Nc-2*Nu)/2, Nsymb); ...     % zero-pad neg. frequencies
              flipud(conj(dataTXm2));...            % data*(-n) (Nu) 
              zeros(1, Nsymb); ...                % 0 at f == 0 (1)
              dataTXm2; ...                         % data(n)  (Nu)
              zeros((Nc-2*Nu)/2 - 1, Nsymb)], 1);   % zero-pad pos. frequencies

% Perform ifft (Nc-IDFT) columnwise
xn = Nc*ifft(Xn, Nc, 1);

% Insert cyclic prefix
xncp = [xn(end-Npre_os+1:end, :); xn]; % insert cyclic prefix

% Parallel to serial
xncp = reshape(xncp, 1, (Nc + Npre_os)*Nsymb); % time-domain ofdm signal w/ cyclic prefix

%% Clipping
% Note: both positive and negative tails are clipped
xncpc = xncp;
xncpc(xncp < 0) = 0; % Clip at zero MUST be done first otherwise distortion will also affect odd subcarriers
xncpc(xncp > sim.rcliptx*sigtx) = sim.rcliptx*sigtx;

% check approximaitons
approx.clip_prob = [qfunc(sim.rcliptx) mean(xncp > sim.rcliptx*sigtx)]; % Clipping probability

%% Quantize at the transmitter
if sim.quantiz
    % With quantization (clip -> quantiz -> interp)  
    % Quantiz
    yqtx = linspace(0, sim.rcliptx*sigtx, 2^sim.ENOB);     % Quantization levels
    delta1 = abs(yqtx(2)-yqtx(1));            % sim.rcliptx*sigtx/(2^sim.ENOB-1);        % (2^sim.ENOB-1) because of clipping done before
    [~, xncpq, varQ1] = quantiz(xncpc, yqtx(1:end-1) + delta1/2, yqtx);

    QuantErrTX = xncpc - xncpq; % used to plot histogram
    
    % Check approximations
    approx.delta1 = [delta1th, delta1];
    approx.varQ1 = [varQ1th, varQ1];
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

%% Ensure signal is positive
% After DA conversion signal might have become negative again. Thus,
% add smal dc bias to make signal positive
if min(xt) < 0
    % Add DC bias
    xtc = xt - min(xt);
else
    disp('No DC bias added!!')
    xtc = xt;
end
% Note: if using Bessel or Gaussian filter this value should be very
% small and thus can be ignored from the penalty calculation.

%% Apply frequency response of the laser
Pt = real(ifft(ifftshift(Hl).*fft(xtc)));  % Note: real() is used to remove residual imag
                                            % part that appears due to numerical error        
                                
% Calculate equivalent transmitted power
Ptx  = mean(Pt);                                     

%% Receiver
% Detect and sample signal 
Its = rx.R*Pt;                      % detect signal

% Add noise
wt = sqrt(rx.Sth*ofdm.fsct)*randn(1, sim.Mct*sim.Nsymb*(Npre_os + Nc));
It = Its + wt;

% Antialiasing filter of the ADC
It = real(ifft(fft(It).*ifftshift(Gadc)));
Its = real(ifft(fft(Its).*ifftshift(Gadc)));
wt = real(ifft(fft(wt).*ifftshift(Gadc)));
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
    % center signal around sigtx/sqrt(2*pi) as defined for quantization.
    % This just removes the additional dc bias that was added at the
    % transmitter to make the signal positive.
    yncp = yncp - mean(yncp) + sigtx/sqrt(2*pi);
    
    % quantization levels
    yqrx = linspace(0, sigtx/sqrt(2*pi) + sim.rcliprx*sigrx, 2^sim.ENOB); % Quantization levels
    delta2 = abs(yqrx(2) - yqrx(1)); 
    
    % Determines what should be the initial point of the quantizer to
    % better use quantizer range.
    p0 = find_opt_p0(yncp, delta2, sigtx/sqrt(2*pi) + sim.rcliprx*sigrx);
    
    yqrx = yqrx + p0*delta2;
        
    % clip before quantization
    yncpc = yncp;
    yncpc(yncp < yqrx(1)) = yqrx(1);
    yncpc(yncp > yqrx(end)) = yqrx(end);
    
    % Quantize
    [~, yncpq, varQ2] = quantiz(yncpc, yqrx(1:end-1) + delta2/2, yqrx);
   
    QuantErrRX = yncpc - yncpq; % used to plot histogram
    
    % Check approximations
    approx.delta2 = [delta2th, delta2];
    approx.varQ2 = [varQ2th, varQ2];
else
    % No quantization
    varQ2 = 0;
    yncpq = yncp;
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

% Demodulate symbols
Yn = fft(yn, Nc, 1)/Nc;             % signal + noise      
Yns = fft(yns, Nc, 1)/Nc;           % signal 
Wn = fft(wn, Nc, 1)/Nc;             % noise

% Calculate SNR at each subcarrier
Yns = mean(abs(Yns).^2, 2);
Yns = Yns(1 + (1:2:Nu)).';      % Only consider positive-frequency subcarriers

Wn = mean(abs(Wn).^2, 2); 
Wn = Wn(1 + (1:2:Nu)).';        % Only consider positive-frequency subcarriers

SNRn = 10*log10(Yns./(Wn + varQ1/Nc*abs(Gch).^2 + varQ2/Nc)); % Measured SNR

% SNRnest is estimated from theory. If BER is estimated using SNRnest, then it
% should be exactly equal to the target BER.
SNRnest = 10*log10(Nc*Pn.*abs(K*Gch).^2./(rx.Sth*ofdm.fs*abs(Gadc(nfc)).^2 + varQ1th*abs(Gch).^2 + varQ2th)); 

% used subcarrier amplitudes (complex conjugate subcarriers are ignored)
dataRXm = Yn(1 + (1:2:Nu), Ndis+1:end-Ndis);            

% AGC caculates what should be the scaling factor to normalize all
% constellations so that dmin = 2 (value expected by qamdemod)
AGCn = sqrt(Pqam)./(K*sqrt(Pn).*Gch);
% Note: Factor of 4 appears because after clipping the amplitude of the
% odd subcarriers (signal bearing subcarriers) are halved.

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
ber.count = sum(numerr)/(sim.Nsymb*sum(bn));

% 95% confidence intervals for the counted BER
[~, interval] = berconfint(numerr, sim.Nsymb*bn);
ber.interval = [min(interval(:,1)), max(interval(:,2))];

% Estimate BER from the SNR at each subcarrier
% Note that to estimate the BER we use SNRn which is calculated from the
% data measured in the simulation. If we had used SNRnest (calculated from
% the theoretical values) the results would match the target BER perfectly
berest = zeros(size(SNRn));
for kk = 1:length(SNRn)
    berest(kk) = berqam(CS(kk), SNRn(kk));
end

ber.est = sum(berest.*bn)/sum(bn);

%% Power
P.Ptx = Ptx;
P.Pn = Pn;
P.Ptxest = Ptxest;

%% Approximations
approx.sigtx = [sigtx std(xncp)];
approx.sigrx = [sigrx std(yncp)];
approx.Ptx = [Ptxest Ptx];

%% Penalties
% Penalty due to nonflat frequency response of the channel (DAC * Hl * ADC)
pp.nonflat = 10*log10(1/sqrt(2*pi)*sqrt(2*sum(K^2./abs(Gch).^2))); 
% Note: the factor K^2 appears multiplying rather than dividing the power
% because the average power depends is given by sigma/(2*pi), where sigma
% is the power before clipping.

% Cyclic prefix
pp.cp = 10*log10(sqrt((Nc+ofdm.Npre_os)/Nc));

% DC-bias
pp.dcbias = 10*log10((mean(xncpc)-min(xt))/1e-3) - 10*log10(mean(xncpc)/1e-3);

%% Show plots if user set verbose to true
if verbose 
    disp('Approximations: [theory, measured]')
    
    approx
   
    plot_stuff
end

% % Right tail of signal ACO-OFDM
% figure
% [n, xout] = hist(yncp, 50);
% dx = abs(xout(1)-xout(2));
% n = n/(trapz(n)*dx);
% ynorm = pdf('normal', xout, sigtx/sqrt(2*pi), sigrx);
% plot(xout, n, xout, ynorm)
% legend('measured', 'approx')
