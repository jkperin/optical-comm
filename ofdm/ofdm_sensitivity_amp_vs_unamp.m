%% Compare the receiver sensitivity of unamplified and amplified IM-DD systems for DC-OFDM
clear, clc

addpath f/
addpath ../f/
addpath ../apd/

sim.OFDM = 'ACO-OFDM'; % either 'DC-OFDM' or 'ACO-OFDM'
sim.Rb = 112e9; % Bit rate
sim.BERtarget = 1.8e-4; % target BER
sim.Mct = 5; % Oversampling ratio to emulate continuous time
Lkm = 0:20; % Fiber length in km

wavelength = 1380e-9; % Transmission wavelength
RIN = -150; % dB/Hz. 
alpha = 0; % chirp parameter of the modulator
modBW = 30e9; % modulator bandwidth
rexdB = -15; % modulator extinction ratio

N0 = (30e-12)^2; % thermal noise PSD

powerAllocation = 'Levin-Campello-MA'; % power allocation method
quantiz = true; % whether quantization is assumed

%% OFDM 
% OFDM constructor: ofdm(Nc, Nu, CS, Rb, power_allocation_type (optional, default = 'palloc')
% Nc : number of subcarriers
% Nu : number of used subcarriers
% CS : nominal constellation size
% Rb : bit rate (b/s)
% power_allocation_type : {'palloc', 'preemphasis'}
if strcmpi(sim.OFDM, 'DC-OFDM')
    ofdm = ofdm(256, 208, 16, 112e9, 'DC', powerAllocation); 
    ENOB = 5;
    rclip = 3; % clipping ratio (controls the DC-bias of the DC-OFDM signal)
elseif strcmpi(sim.OFDM, 'ACO-OFDM')
    ofdm = ofdm(256, 208, 64, 112e9, 'ACO', powerAllocation); 
    ENOB = 6;
    rclip = 3.3;
end

ofdm.set_cyclic_prefix(5, 5); % set cyclic prefix length. Should be consistent with channel memory length
noiseBW = ofdm.fs/2;
varRIN =   @(Prx) 10^(RIN/10)*Prx.^2*noiseBW;

sim.fs = ofdm.fs*sim.Mct;

% Modulator
filt = design_filter('butter', 5, noiseBW/(sim.fs/2));
Nhold = sim.Mct;
hZOH = 1/Nhold*ones(1, Nhold);
Hdac = filt.H(ofdm.fc/sim.fs).*freqz(hZOH, 1, ofdm.fc, sim.fs).*exp(1j*2*pi*ofdm.fc/sim.fs*(Nhold-1)/2);

Mod.filt = design_filter('two-pole', modBW, sim.fs);
% Tx.Mod.filt = design_filter('butter', 5, Tx.Mod.BW/(sim.fs/2));
Hmod = Mod.filt.H(ofdm.fc/sim.fs);

% modfc = modBW/sqrt(sqrt(2)-1); % converts to relaxation frequency
% Hmod = 1./(1 + 2*1j*ofdm.fc/modfc - (ofdm.fc/modfc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)

% Fiber
Fiber = fiber(0);

% Optical amplifier
OptAmp = OpticalAmplifier('ConstantGain', 20, 5, wavelength);
% OptAmp = OpticalAmplifier('ConstantOutputPower', 0, 5, wavelength);
BWopt = sim.fs;
Npol =2;

% Photodiode
PD = pin(1, 10e-9);
Hpd = PD.H(ofdm.fc);

% APD
% constructor: apd(GaindB, ka, BW, R, Id)
APD = apd(10, 0.18, [24e9 290e9], 0.74, 40e-9);
% APD = apd(10, 0.18, Inf, 0.74, 40e-9);
Hapd = APD.H(ofdm.fc);

% ADC
Hadc = filt.H(ofdm.fc/sim.fs);

%% Build ananonymous function for the quantization noise generated at the transmitter and receiver
% Quantization
if quantiz
    varQuantizTx = ofdm.quantization_noise_var(ENOB, rclip); % quantization at the transmitter
    varQuantizRx = ofdm.quantization_noise_var(ENOB, rclip); % quantization at the transmitter
else
    varQuantizTx = @(Psig) 0;
    varQuantizRx = @(Psig) 0;
end

PrecUnamp = zeros(size(Lkm));
PrecAPDopt = zeros(size(Lkm));
PrecAmp = zeros(size(Lkm));
OSNRdB = zeros(size(Lkm));
D = zeros(size(Lkm));
maxTol = 1e-3; % max tolerance for power convergence
maxIterations = 50; % max number of iterations
warning('on');
for k = 1:length(Lkm)
    fprintf('L = %d\n', Lkm(k));
    Fiber.L = 1e3*Lkm(k);
    D(k) = 1e6*Fiber.D(wavelength)*Lkm(k);
       
    % Generate new copies of OFDM class for un- and amplified simulations
    [ofdmUnamp, ofdmAPD, ofdmAmp] = ofdm.copy();   
    
    %% Amplified system
    % Channel frequency response
    Hch = Hdac.*Hmod.*Fiber.Himdd(ofdm.fc, wavelength, alpha, 'small signal').*Hpd.*Hadc;
    
    tol = Inf; % tolerance
    Prx = 1e-6; % starting received power
    n = 0; % number of iterations
    while tol > maxTol && n < maxIterations
        % Note: The power_allocation method of the class OFDM assumes that
        % the noise variance does not depend on the power. However, in this
        % case the noise is signal dependent. Hence, an iterative process
        % must be performed in order to estimate the correct power: 
        % 1. The received power is given an certain value (Prx = 1e-3).
        % 2. The noise variance is calculated assuming that value of Prx.
        % 3. power_allocation method is called for that noise variance. The
        % power allocation algorithm calcultes how much power must be
        % allocated in each subcarrier to achieve the target BER. The
        % power in all subcarriers determine the average power, which is
        % calculated by the method dc_bias of the class OFDM.
        % 4. If the new power calculated in the last step differs from the
        % Prx used to calculate the noise variance in step 2, then this
        % procedure is repeated until the allocated power (step 3) agrees 
        % with power used to calculate the noise variance (step 2) to a
        % predefined tolerance.
        
        % Note: Ppd is the signal power at the input of the photodiode, 
        % while Prx is the signal power at the input of the optical amplifier
        if strcmpi(OptAmp.Operation, 'ConstantOutputPower')
            Ppd = dBm2Watt(OptAmp.outputPower); % power at the photodiode (after the amplifier)
        else
            Ppd = Prx*OptAmp.Gain;
        end

        varNoiseAmp = 1/ofdm.Nc*(abs(Hadc).^2*noiseBW*N0... % thermal 
                   + abs(Hadc.*Hpd).^2*PD.varShot(Ppd, noiseBW)... % shot
                   + abs(Hch*PD.R).^2*varRIN(Ppd)... % intensity noise
                   + abs(Hpd.*Hadc*PD.R).^2*(OptAmp.varNoiseDD(Prx, noiseBW, BWopt, Npol))... % sig-spont + spont-spont
                   + abs(Hch).^2*varQuantizTx(PD.R*Ppd) + varQuantizRx(PD.R*Ppd)); % quantization noise 
                
        ofdmAmp.power_allocation(Hch, varNoiseAmp, sim.BERtarget, ~true);
        [~, Ppdnew] = ofdmAmp.dc_bias(ofdmAmp.Pn, rclip, Hdac, rexdB); % Power required at receiver DSP i.e., after amplification
        Prxnew = Ppdnew/(PD.Geff*OptAmp.Gain); % Refer power to the input of amplifier input
        tol = abs(Prx - Prxnew)/Prxnew;
        Prx = Prxnew;
        n = n + 1;
    end
    if n == maxIterations
        PrecAmp(k) = NaN;
        OSNRdB(k) = NaN;
        disp('Receiver sensitivity calculation did not converge for amplified system');
    else
        PrecAmp(k) = Watt2dBm(Prx);  % receiver sensitivity
        OSNRdB(k) = 10*log10(Ppd/(2*OptAmp.Ssp*OptAmp.BWref));
        Varamp = PD.R^2*(OptAmp.varNoiseDD(Prx, noiseBW, BWopt, Npol));
        disp('-- Amplified')
        fprintf('Thermal/Sig-ASE = %f dB\nShot/Sig-ASE = %f dB\nRIN/Sig-ASE =%f dB\nQuant/Sig-ASE = %f dB\n',...
            10*log10(noiseBW*N0/Varamp), 10*log10(PD.varShot(Ppd, noiseBW)/Varamp), 10*log10(varRIN(Ppd)/Varamp), 10*log10((varQuantizTx(Ppd) + varQuantizRx(Ppd))/Varamp))
        
        % Approximation
%         Pn = ofdmAmp.Pn; % subcarrier powers referred to receiver output
%         [Padj, ~, exitflag] = fzero(@(Padj) log10(ofdmAmp.calc_ber(10*log10(ofdm.Nc*Padj*Pn...
%             /(ofdm.fs*(2*OptAmp.Gain*PD.R*Prx*OptAmp.Ssp)...
%             + 2*rclip^2*ofdmAmp.var/(3*2^(2*ENOB))...
%             + 0*varRIN(OptAmp.Gain*PD.R*Prx)))) - log10(sim.BERtarget)), 1);
%         PrecAmpApprox(k) = 2*sum(Padj*Pn/(PD.Geff*OptAmp.Gain)^2);
%         1;
    end
    
    %% APD-based system: shot and thermal noise comparable
    % Gain optimization
    GainAPD = 1:15;
    PrecAPDGain = zeros(size(GainAPD));
    for kk = 1:length(GainAPD)
        APD.Gain = GainAPD(kk);
        Hapd = APD.H(ofdm.fc); % photodiode frequency response
        Hch = Hdac.*Hmod.*Fiber.Himdd(ofdm.fc, wavelength, alpha, 'small signal').*Hapd.*Hadc;  % Channel frequency response
        
        tol = Inf; % tolerance
        Prx = 1e-6; % starting received power
        n = 0; % number of iterations
        while tol > maxTol && n < maxIterations
            % Note: the same iterative method described for amplified systems
            % is applied here.

            varNoiseAPD = 1/ofdm.Nc*(abs(Hadc).^2*noiseBW*N0... % thermal noise 
                    + abs(Hadc.*Hpd).^2*APD.varShot(Prx, noiseBW)... % shot
                    + abs(Hch).^2*varRIN(APD.Geff*Prx)... % intensity noise
                    + abs(Hch).^2*varQuantizTx(APD.Geff*Prx) + varQuantizRx(APD.Geff*Prx)); % quantization noise           

            ofdmAPD.power_allocation(Hch, varNoiseAPD, sim.BERtarget, ~true);
            [~, Prxnew] = ofdmAPD.dc_bias(ofdmAPD.Pn, rclip, Hdac, rexdB);
            Prxnew = Prxnew/APD.Geff; % Refer power to the input of the APD
            tol = abs(Prx - Prxnew)/Prxnew;
            Prx = Prxnew;
            n = n + 1;
        end
        if n == maxIterations
            PrecAPDGain(kk) = NaN;
            fprintf('Receiver sensitivity calculation did not converge for APD-based system with APD gain = %d\n', APD.Gain);
        else
            PrecAPDGain(kk) = Prx;
        end
    end
    
    % Calculate performance at optimal gain
    PrecAPDGain = Watt2dBm(PrecAPDGain);
    g = linspace(GainAPD(1), GainAPD(end));
    Prec = interp1(GainAPD, PrecAPDGain, g);    
    [PrecAPDopt(k), idx] = min(Prec);
    APD.Gain = g(idx);
    
    figure(100), hold on, box on
    plot(g, Prec, 'LineWidth', 2)
    plot(APD.Gain, PrecAPDopt(k), '*r')
    xlabel('APD gain')
    ylabel('Receiver sensitivity (dBm)')
    drawnow
    
    Prec = dBm2Watt(PrecAPDopt(k));
    varShot = PD.varShot(Prec, noiseBW);
    disp('-- APD')
    fprintf('Optimal gain = %.2f\n', APD.Gain)
    fprintf('Thermal/Shot = %f dB\nRIN/Shot = %f dB\nQuantiz/Shot = %f dB\n',...
    10*log10(N0*noiseBW/varShot), 10*log10(varRIN(APD.Geff*Prec)/varShot),...
    10*log10((varQuantizTx(APD.Geff*Prec) + varQuantizRx(APD.Geff*Prec))/varShot))
    
    %% Unamplified system: thermal-noise dominant
    % Channel frequency response
    Hch = Hdac.*Hmod.*Fiber.Himdd(ofdm.fc, wavelength, alpha, 'small signal').*Hpd.*Hadc;
    
    tol = Inf; % tolerance
    Prx = 1e-3; % starting received power
    n = 0; % number of iterations
    while tol > maxTol && n < maxIterations
        % Note: the same iterative method described for amplified systems
        % is applied here.
        
        varNoiseUnamp = 1/ofdm.Nc*(abs(Hadc).^2*noiseBW*N0... % thermal noise 
                + abs(Hadc.*Hpd).^2*PD.varShot(Prx, noiseBW)... % shot
                + abs(Hch).^2*varRIN(PD.Geff*Prx)... % intensity noise
                + abs(Hch).^2*varQuantizTx(PD.Geff*Prx) + varQuantizRx(PD.Geff*Prx)); % quantization noise        
    
        ofdmUnamp.power_allocation(Hch, varNoiseUnamp, sim.BERtarget, ~true);
        [~, Prxnew] = ofdmUnamp.dc_bias(ofdmUnamp.Pn, rclip, Hdac, rexdB);
        Prxnew = Prxnew/PD.Geff; % Refer power to the input of the PD
        tol = abs(Prx - Prxnew)/Prxnew;
        Prx = Prxnew;
        n = n + 1;
    end
    if n == maxIterations
        PrecUnamp(k) = NaN;
        disp('Receiver sensitivity calculation did not converge for unamplified');
    else
        PrecUnamp(k) = Watt2dBm(Prx);
        Vartherm = noiseBW*N0;
        disp('-- Unamplified')
        fprintf('Shot/Thermal = %f dB\nRIN/Thermal = %f dB\nQuantiz/Thermal = %f dB\n',...
        10*log10(PD.varShot(Prx, noiseBW)/Vartherm), 10*log10(varRIN(Prx)/Vartherm), 10*log10((varQuantizTx(Prx) + varQuantizRx(Prx))/Vartherm))
    end
end

figure(1), hold on, box on
hplot(1) = plot(D, PrecUnamp, '-', 'LineWidth', 2, 'DisplayName', [sim.OFDM ' unamplified']);
hplot(2) = plot(D, PrecAPDopt, '--', 'Color', get(hplot(1), 'Color'), 'LineWidth', 2, 'DisplayName', [sim.OFDM ' APD']);
hplot(3) = plot(D, PrecAmp, ':', 'Color', get(hplot(1), 'Color'), 'LineWidth', 2, 'DisplayName', [sim.OFDM ' amplified']);
% plot(D, PrecAmpApprox, '-', 'Color', get(hplot(3), 'Color'), 'LineWidth', 2, 'DisplayName', [sim.OFDM ' amplified approx']);
legend('-dynamiclegend')
xlabel('Dispersion (ps/nm)', 'FontSize', 12)
ylabel('Receiver sensitivity (dBm)', 'FontSize', 12)
set(gca, 'FontSize', 12)

figure(2), hold on, box on
plot(D, OSNRdB, '-', 'Color', get(hplot(3), 'Color'), 'LineWidth', 2, 'DisplayName', sim.OFDM)
legend('-dynamiclegend')
xlabel('Dispersion (ps/nm)', 'FontSize', 12)
ylabel('OSNR required (dB)', 'FontSize', 12)
set(gca, 'FontSize', 12)

% Receiver sensitivity for an ideal ISI-free IM-DD channel
M = ofdm.CS;
BERaco = @(P) 4/log2(M)*(1 - 1/sqrt(M))*qfunc(P/rclip*sqrt((3/(M-1))*log2(M)/(N0*sim.Rb)));
PrxdBm = fzero(@(P) log10(BERaco(dBm2Watt(P))) - log10(sim.BERtarget), -20)
