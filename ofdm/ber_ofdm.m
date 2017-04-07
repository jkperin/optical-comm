function [ber, ofdm, OSNRdB] = ber_ofdm(ofdm, Tx, Fibers, Rx, sim)
%% Calculate BER of DC or ACO-OFDM in IM-DD system
% This function calculates appropriate cyclic prefix length given the
% channel memory length. It calls ber_dc_ofdm_montecarlo, which performs
% Montecarlo simulation to obtain the BER at each transmitted power levels.
% ber_dc_ofdm_montecarlo also estimates theoretical BER based on Gaussian
% approximation.
% Inputs:
% - ofdm: OFDM class
% - Tx: struct with transmitter paramters
% - Fibers: array of Fiber classes containing the fibers used in
% transmittion i.e., Fibers = [SMF DCF] or Fibers = [SMF];
% - Rx: struct with receiver parameters
% - sim: struct with simulation parameters
% Outputs:
% - ber: BER using Montecarlo simulation and Gaussian approximation
% - ofdm: class ofdm after power allocation
% - OSNRdB: OSNR. Only meaninful if simulation includes optical amplifier

%% Calculate components frequency response in order to estimate cyclic prefix length
% Total (residual) dispersion
Dtotal = 0;
linkAtt = 1;
for k = 1:length(Fibers)
    fiberk = Fibers(k);
    Dtotal = Dtotal + fiberk.D(Tx.Laser.wavelength)*fiberk.L;
    linkAtt = linkAtt*fiberk.link_attenuation(Tx.Laser.wavelength);
end
fprintf('Total dispersion at %.2f nm: %.3f ps/nm\n', Tx.Laser.lambda*1e9, 1e3*Dtotal)

% Create fiber object with equivalent dispersion.
equivFiber = fiber(1, @(lamb) 0, @(lamb) Dtotal);
Hfiber = @(f) equivFiber.Himdd(f, Tx.Laser.wavelength, Tx.alpha, 'small signal'); % frequency response of the channel (used in designing the equalizer)

% ZOH and DAC frequency responses
Nhold = sim.Mct/Tx.DAC.ros;
hZOH = 1/Nhold*ones(1, Nhold);
Hdac = @(f, fs) Tx.DAC.filt.H(f/fs).*freqz(hZOH, 1, f, fs)...
    .*exp(1j*2*pi*f/fs*(Nhold-1)/2);

% Modulator frequency response
if isfield(Tx, 'Mod')
    Hmod = @(f, fs) Tx.Mod.filt.H(f/fs); % group delay was already removed
else
    Hmod = @(f, fs) ones(size(f));
end

% PD
Hpd = @(f) Rx.PD.H(f);

% Antialiasing filter
Hadc = @(f, fs) Rx.ADC.filt.H(f/fs);

%% Calculate cyclic prefix
Hch_fun = @(f, fs) Hdac(f, fs).*Hmod(f, fs).*Hfiber(f).*Hpd(f).*Hadc(f, fs);

if isempty(ofdm.Nneg) % calculate cyclic prefix if not yet calculated
    ofdm.cyclic_prefix(@(f, fs) Hch_fun, sim.shouldPlot('Cyclic prefix'));
    fprintf('> Cyclic prefix:\nNneg = %d\nNpos = %d\n', ofdm.Nneg, ofdm.Npos);
end

Tx.Hdac = Hdac(ofdm.fc, sim.fs); % used in DC-bias calculation
Hmod = Hmod(ofdm.fc, sim.fs); % modulator
Hadc = Hadc(ofdm.fc, sim.fs); 
Hpd = Hpd(ofdm.fc); % photodiode's frequency response at the ofdm subcarriers frequency
Hch = Hch_fun(ofdm.fc, sim.fs);
Hch_pd = Tx.Hdac.*Hmod.*Hfiber(ofdm.fc).*Hadc; % channel response without PD

% Noise bandwidth
noiseBW = ofdm.fs/2;

% Thermal noise
varTherm = Rx.N0*noiseBW; % variance of thermal noise
varRIN = @(Prx) 0;
if isfield(sim, 'RIN') && sim.RIN
    varRIN =   @(Prx) 10^(Tx.RIN/10)*Prx.^2*noiseBW;
end

% Quantization
if isfield(sim, 'quantiz') && sim.quantiz
    varQuantizTx = ofdm.quantization_noise_var(Tx.DAC.resolution, Tx.rclip); % quantization at the transmitter
    varQuantizRx = ofdm.quantization_noise_var(Rx.ADC.ENOB, Tx.rclip); % quantization at the transmitter
else
    varQuantizTx = @(Psig) 0;
    varQuantizRx = @(Psig) 0;
end

%% Calculate BER
Ptx = dBm2Watt(Tx.PtxdBm); % Transmitted power

ber.gauss = zeros(size(Ptx));
ber.count = zeros(size(Ptx));
ber.theory = zeros(size(Ptx));
OSNRdB = zeros(size(Ptx));
for k = 1:length(Ptx)
    Tx.Ptx = Ptx(k);
           
    % Estimate noise to perform power allocation
    if isfield(sim, 'preAmp') && sim.preAmp % amplified system: signal-spontaneous beat noise dominant
        % Prx = signal power referred to the photodiode's input
        if strcmpi(Rx.OptAmp.Operation, 'ConstantOutputPower')
            Prx = dBm2Watt(Rx.OptAmp.outputPower); % Received power is constant at the amplifier output
            Rx.OptAmp.Gain = Prx/(linkAtt*Tx.Ptx);
        else
            Prx = linkAtt*Rx.OptAmp.Gain*Tx.Ptx;
        end
        BWopt = sim.fs; % optical filter (OF) noise bandwidth = sampling rate, since OF is not included 
        % Not divided by 2 because optical filter is a bandpass filter

        Npol = 2; % number of polarizations. Npol = 1, if polarizer is present, Npol = 2 otherwise.
        varNoise = 1/ofdm.Nc*(abs(Hadc).^2*varTherm... % thermal 
                           + abs(Hadc.*Hpd).^2*Rx.PD.varShot(Prx, noiseBW)... % shot
                           + abs(Rx.PD.R*Hch).^2*varRIN(Prx)... % intensity noise
                           + abs(Rx.PD.R*Hpd.*Hadc).^2*(Rx.OptAmp.varNoiseDD(Prx/Rx.OptAmp.Gain, noiseBW, BWopt, Npol))... % sig-spont + spont-spont
                           + abs(Hch).^2*varQuantizTx(Rx.PD.R*Prx) + abs(Hadc).^2*varQuantizRx(Rx.PD.R*Prx)); % quantization noise 
        % Note: Prx is divided by amplifier gain to obtain power at the amplifier input
    else % unamplified system: photodiode is either PIN or APD
        Prx = Tx.Ptx*linkAtt; % power referred to the photodiode's input
        
        % if PD == APD, perform gain optimization
        if strcmpi(class(Rx.PD), 'APD') && isfield(sim, 'optimizeAPDGain') && sim.optimizeAPDGain
            disp('Optimizing APD gain...')
            [optGain, ~, exitflag] = fminbnd(@(Gain) iterate(Gain, Tx, ofdm, Rx, sim), 1, 30); % limit gain search between 1 and 30
            
            if exitflag ~= 1
                warning('APD gain optimization exited with exitflag = %d', exitflag)
            end
            fprintf('Optimal APD gain = %.2f\n', optGain);
            Rx.PD.Gain = optGain;
            Hch = Hch_pd.*Rx.PD.H(ofdm.fc);
        end
        
        varNoise = 1/ofdm.Nc*(abs(Hadc).^2*varTherm... % thermal noise 
                            + abs(Hadc.*Hpd).^2*Rx.PD.varShot(Prx, noiseBW)... % shot
                            + abs(Hch).^2*varRIN(Rx.PD.Geff*Prx)... % intensity noise
                            + abs(Hch).^2*varQuantizTx(Rx.PD.Geff*Prx) + abs(Hadc).^2*varQuantizRx(Rx.PD.Geff*Prx)); % quantization noise                   
    end    
    
    % Power allocation
    ofdm.power_allocation(Hch, varNoise, sim.BERtarget, sim.shouldPlot('Power allocation'));
    % Note: result of power allocation is referred to the receiver, after
    % amplification, photodetection, etc. It is not necessary to refer the
    % power to the transmitter because only the relative powers of the
    % subcarriers is important. The total power is limited by the
    % transmitted power
    
    % Check if power is at reasonable levels
    assert(sqrt(ofdm.var) < 1, 'ber_ofdm: Power is too high: %.2G (mW)\nCannot improve performance by increasing power.\n', 1e3*sqrt(2*sum(ofdm.Pn)))
    
    % Theoretical BER
    Padj = ofdm.adjust_power_allocation(ofdm.Pn.*abs(Tx.Hdac).^2, Prx*Rx.PD.Geff, Tx.rclip, Tx.rexdB); % Factor to refer subcarriers powers to the receiver
    ber.theory(k) = ofdm.calc_ber(10*log10(abs(ofdm.K*Hch).^2.*ofdm.Pn*Padj./varNoise));
    
    % Montecarlo simulation
    sim.Hch = Hch; % frequency response of the channel at each subcarrier
    sim.varNoise = varNoise; % noise variance in each subcarrier
    [ber.count(k), ber.gauss(k), ~, OSNRdB(k)] = ber_ofdm_montecarlo(ofdm, Tx, Fibers, Rx, sim);
    
    fprintf('Ptx = %.2f dBm\n- BER(theory) = %g\n- BER(Gauss) = %g\n- BER(count) = %g\n',...
        Tx.PtxdBm(k), ber.theory(k), ber.gauss(k), ber.count(k));
    
    % Whether to stop simulation when counted BER reaches 0
    if sim.stopSimWhenBERReaches0 && ber.count(k) == 0
        break;
    end 
end

%% Plots
if sim.shouldPlot('BER') && length(ber.count) > 1
    figure(1), hold on, box on
    hline = plot(Tx.PtxdBm, log10(ber.count), '-o', 'LineWidth', 2);
    plot(Tx.PtxdBm, log10(ber.gauss), '-', 'Color', get(hline, 'Color'), 'LineWidth', 2)
    plot(Tx.PtxdBm, log10(ber.theory), ':k', 'LineWidth', 2)
    legend('Counted', 'Gaussian approx.', 'Theory')
    axis([Tx.PtxdBm(1) Tx.PtxdBm(end) -8 0])
    xlabel('Received power (dBm)', 'FontSize', 12)
    ylabel('log_{10}(BER)', 'FontSize', 12)
    set(gca, 'FontSize', 12)
    drawnow
end

function Psens = iterate(Gain, Tx, ofdm, Rx, sim)
%% Function used in APD gain optimization
% For each value of gain, calculate new noise variance, do power
% allocaiton, and determine new average optical power
    Rx.PD.Gain = Gain;
    Hpd = Rx.PD.H(ofdm.fc); % photodiode's frequency response at the ofdm subcarriers frequency
    Hch_new = Hch_pd.*Hpd;
    
    tol = Inf; % tolerance
    Prec = 1e-6; % starting received power
    n = 0; % number of iterations
    maxTol = 1e-3; % max tolerance for power convergence
    maxIterations = 50; % max number of iterations
    while tol > maxTol && n < maxIterations
        % Note: the same iterative method described for amplified systems
        % is applied here.

        varNoiseAPD = 1/ofdm.Nc*(abs(Hadc).^2*varTherm... % thermal noise 
                               + abs(Hpd.*Hadc).^2*Rx.PD.varShot(Prec, noiseBW)... % shot
                               + abs(Hch_new).^2*varRIN(Rx.PD.Geff*Prec)... % intensity noise
                               + abs(Hch_new).^2*varQuantizTx(Rx.PD.Geff*Prec) + abs(Hadc).^2*varQuantizRx(Rx.PD.Geff*Prec)); % quantization noise

        ofdm.power_allocation(Hch_new, varNoiseAPD, sim.BERtarget, false);
        [~, Prxnew] = ofdm.dc_bias(ofdm.Pn, Tx.rclip, Tx.Hdac, Tx.rexdB);
        Prxnew = Prxnew/Rx.PD.Geff;
        tol = abs(Prec - Prxnew)/Prxnew;
        Prec = Prxnew;
        n = n + 1;
    end
    
    Psens = Watt2dBm(Prec);
end
end