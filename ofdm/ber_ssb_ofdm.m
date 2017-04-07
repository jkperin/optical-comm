function [ber, ofdm, OSNRdB] = ber_ssb_ofdm(ofdm, Tx, Fibers, Rx, sim)
%% Calculate BER of SSB-OFDM in IM-DD system
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
CSPR = 10^(Tx.CSPRdB/10);

ber.gauss = zeros(size(Ptx));
ber.count = zeros(size(Ptx));
ber.theory = zeros(size(Ptx));
OSNRdB = zeros(size(Ptx));
for k = 1:length(Ptx)
    Tx.Ptx = Ptx(k);
           
    % Estimate noise to perform power allocation
    % Prx = signal power referred to the photodiode's input
    if strcmpi(Rx.OptAmp.Operation, 'ConstantOutputPower')
        Prx = dBm2Watt(Rx.OptAmp.outputPower); % Received power is constant at the amplifier output
        Rx.OptAmp.Gain = Prx/(linkAtt*Tx.Ptx);
    else
        Prx = linkAtt*Rx.OptAmp.Gain*Tx.Ptx;
    end
    
    Ps = (Prx/Rx.OptAmp.Gain)/(1 + CSPR);
    Pc = CSPR*Ps;
    
    varNoise = 1/(ofdm.Nc)*(abs(Hadc).^2*varTherm... % thermal 
                       + abs(Hadc.*Hpd).^2*Rx.PD.varShot(Prx, noiseBW)... % shot
                       + abs(Rx.PD.R*Hch).^2*varRIN(Prx)... % intensity noise
                       + abs(Rx.PD.R*Hpd.*Hadc).^2*(4*Rx.OptAmp.Gain*(Ps + Pc)*Rx.OptAmp.Ssp*noiseBW)... % sig-spont + spont-spont
                       + varQuantizRx(Rx.PD.R*Rx.OptAmp.Gain*sqrt(2*Pc*Ps))); % quantization noise 
    % Note: Prx is divided by amplifier gain to obtain power at the amplifier input 
    
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
    Pnnorm = ofdm.Pn*Ps/sum(ofdm.Pn); % adjust for desired power 
    ber.theory(k) = ofdm.calc_ber(10*log10(abs(Hch).^2.*Pnnorm./varNoise));
    
    % Montecarlo simulation
    sim.Hch = Hch; % frequency response of the channel at each subcarrier
    sim.varNoise = varNoise; % noise variance in each subcarrier
    [ber.count(k), ber.gauss(k), ~, OSNRdB(k)] = ber_ssb_ofdm_montecarlo(ofdm, Tx, Fibers, Rx, sim);
    
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
