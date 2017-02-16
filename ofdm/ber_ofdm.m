function [ber, ofdm, SNRdB, OSNRdB] = ber_ofdm(ofdm, Tx, Fibers, Rx, sim)
%% Calculate BER of IM-DD system
% This function calculates appropriate cyclic prefix length given the
% channel memory length. It calls ber_dc_ofdm_montecarlo, which performs
% Montecarlo simulation to obtain the BER at each transmitted power levels.
% ber_dc_ofdm_montecarlo also estimates theoretical BER based on Gaussian
% approximation.
% Inputs:
% - mpam: OFDM class
% - Tx: struct with transmitter paramters
% - Fibers: array of Fiber classes containing the fibers used in
% transmittion i.e., Fibers = [SMF DCF] or Fibers = [SMF];
% - Rx: struct with receiver parameters
% - sim: struct with simulation parameters
% Outputs:
% - ber: BER using Montecarlo simulation and Gaussian approximation
% - ofdm: class ofdm after power allocation
% - SNRdB: average SNRdB
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
Hfiber = @(f) equivFiber.Himdd(f, Tx.Laser.wavelength, Tx.alpha, 'large signal'); % frequency response of the channel (used in designing the equalizer)

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

% Antialiasing filter
Hadc = @(f, fs) Rx.ADC.filt.H(f/fs);

%% Calculate cyclic prefix
Hch = @(f, fs) Hdac(f, fs).*Hmod(f, fs).*Hfiber(f).*Hadc(f, fs);

if isempty(ofdm.Nneg) % calculate cyclic prefix if not yet calculated
    ofdm.cyclic_prefix(Hch, sim.shouldPlot('Cyclic prefix'));
    fprintf('> Cyclic prefix:\nNneg = %d\nNpos = %d\n', ofdm.Nneg, ofdm.Npos);
end

Tx.Hdac = Hdac(ofdm.fc, sim.fs); % used in DC-bias calculation
Hadc = Hadc(ofdm.fc, sim.fs); 
Hch = Hch(ofdm.fc, sim.fs); % channel frequency response at the subcarriers frequency

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
SNRdB = zeros(size(Ptx));
OSNRdB = zeros(size(Ptx));
for k = 1:length(Ptx)
    Tx.Ptx = Ptx(k);
           
    % Estimate noise to perform power allocation
    if isfield(sim, 'preAmp') && sim.preAmp % amplified system: signal-spontaneous beat noise dominant
        Prx = dBm2Watt(Rx.OptAmp.outputPower); % Received power is constant at the amplifier output
        BWopt = sim.fs; % optical filter (OF) noise bandwidth = sampling rate, since OF is not included 
        % Not divided by 2 because optical filter is a bandpass filter

        Npol = 2; % number of polarizations. Npol = 1, if polarizer is present, Npol = 2 otherwise.
        varNoise = 1/ofdm.Nc*(abs(Hadc).^2)*(varTherm + Rx.PD.varShot(Prx, noiseBW) + varRIN(Prx)... % thermal + shot + RIN
                    + Rx.PD.R^2*(Rx.OptAmp.varNoiseDD(Tx.Ptx, noiseBW, BWopt, Npol)))... % sig-spont + spont-spont
                        + 1/ofdm.Nc*(abs(Hch).^2*varQuantizTx(Prx) + varQuantizRx(Prx)); % quantization noise 
        % Note: Prx is divided by amplifier gain to obtain power at the amplifier input
    else % unamplified system: thermal-noise dominant
        Prx = Tx.Ptx*linkAtt;
        
        varNoise = 1/ofdm.Nc*(abs(Hadc).^2)*(varTherm... % thermal noise 
            + Rx.PD.varShot(Prx, noiseBW) + varRIN(Prx))... % shot noise + RIN
            + 1/ofdm.Nc*(abs(Hch).^2*varQuantizTx(Prx) + varQuantizRx(Prx)); % quantization noise
    end    
    
    % Power allocation
    ofdm.power_allocation(Hch, varNoise, sim.BERtarget, sim.shouldPlot('Power allocation'));
    % Check if power is at reasonable levels
    assert(sqrt(ofdm.var) < 1, 'ber_ofdm: Power is too high: %.2G (mW)\nCannot improve performance by increasing power.\n', 1e3*sqrt(2*sum(ofdm.Pn)))
    
    % Theoretical BER
    Padj = ofdm.adjust_power_allocation(ofdm.Pn.*abs(Tx.Hdac).^2, Prx, Tx.rclip, Tx.rexdB);
    ber.theory(k) = ofdm.calc_ber(10*log10(abs(ofdm.K*Hch).^2.*ofdm.Pn*Padj./varNoise));
    % Note: The scaling factor (Prx/Pmean) is squared because the average
    % power is proportional to sqrt(ofdm.var) = sqrt(2*sum(Pn))
    
    % Montecarlo simulation
    sim.Hch = Hch; % frequency response of the channel at each subcarrier
    sim.varNoise = varNoise; % noise variance in each subcarrier
    [ber.count(k), ber.gauss(k), SNRndB, OSNRdB(k)] = ber_ofdm_montecarlo(ofdm, Tx, Fibers, Rx, sim);
    SNRdB(k) = sum(SNRndB.*log2(ofdm.CS))/sum(log2(ofdm.CSn));
    
    fprintf('Ptx = %.2f dBm\n- BER(theory) = %g\n- BER(Gauss) = %g\n- BER(count) = %g\n',...
        Tx.PtxdBm(k), ber.theory(k), ber.gauss(k), ber.count(k));
    
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