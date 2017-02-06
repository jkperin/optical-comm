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

% ZOH and DAC frequency response
Nhold = sim.Mct/Tx.DAC.ros;
hZOH = 1/Nhold*ones(1, Nhold);
Hdac = @(f, fs) Tx.DAC.filt.H(f/fs).*freqz(hZOH, 1, f/fs)...
    .*exp(1j*2*pi*f/fs*(Nhold-1)/2);

% Modulator frequency response
if isfield(Tx, 'Mod')
    Hmod = @(f, fs) Tx.Mod.filt.H(f/fs); % group delay was already removed
else
    Hmod = @(f, fs) ones(size(f));
end

% Antialiasing filter
Hrx = @(f, fs) Rx.ADC.filt.H(f/fs);

%% Calculate cyclic prefix
Hch = @(f, fs) Hdac(f, fs).*Hmod(f, fs).*Hfiber(f).*Hrx(f, fs);

if isempty(ofdm.Nneg) % calculate cyclic prefix if not yet calculated
    ofdm.cyclic_prefix(Hch, sim.shouldPlot('Cyclic prefix'));
    fprintf('> Cyclic prefix:\nNneg = %d\nNpos = %d\n', ofdm.Nneg, ofdm.Npos);
end

Tx.Mod.Hdac = Hdac(ofdm.fc, sim.fs); % used in DC-bias calculation
Hch = Hch(ofdm.fc, sim.fs); % channel frequency response at the subcarriers frequency

% Noise bandwidth
noiseBW = ofdm.fs/2;

% Thermal noise
varTherm = Rx.N0*noiseBW; % variance of thermal noise

% Check if power is at reasonable levels
assert(sqrt(ofdm.var) < 1e-3, 'Power is too high: %.2G (mW)\nCannot improve performance by increasing power.\n', 1e3*sqrt(2*sum(ofdm.Pn)))

% Quantization

if isfield(sim, 'quantiz') && sim.quantiz
    % At the receiver Prx = sigmatx*rclip, where sigmatx = sqrt(2*sum(Pn));
    % This ensures that quantization noise is normalized to the
    % received power, which is what the noise calculations assume.
    varQuantizTx = @(Psig) (1 - 2*qfunc(Tx.rclip))*(2*Psig/(2^Tx.DAC.resolution-1))^2/12; % quantization at the transmitter
    varQuantizRx = @(Psig) (1 - 2*qfunc(Tx.rclip))*(2*Psig/(2^Rx.ADC.ENOB-1))^2/12; % quantization at the transmitter
    % Note: this assumes that quantization noise from the transmitter
    % has the same effect as the quantization noise from the receiver.
    % This is not true, however, since the quantization noise from the
    % transmitter is also filtered by the channel.
else
    varQuantizTx = @(Psig) 0;
    varQuantizRx = @(Psig) 0;
end
    

%% Calculate BER
Ptx = dBm2Watt(Tx.PtxdBm); % Transmitted power

ber.gauss = zeros(size(Ptx));
ber.count = zeros(size(Ptx));
SNRdB = zeros(size(Ptx));
OSNRdB = zeros(size(Ptx));
for k = 1:length(Ptx)
    Tx.Ptx = Ptx(k);
           
    % Estimate noise to perform power allocation
    if isfield(sim, 'preAmp') && sim.preAmp % amplified system: signal-spontaneous beat noise dominant
        Prx = dBm2Watt(Rx.OptAmpOutPowerdBm);
        BWopt = sim.fs; % optical filter (OF) noise bandwidth = sampling rate, since OF is not included 
        % Not divided by 2 because optical filter is a bandpass filter

        % Noise std for intensity level Plevel
        Npol = 2; % number of polarizations. Npol = 1, if polarizer is present, Npol = 2 otherwise.
        varNoise = 1/ofdm.Nc*ones(1, ofdm.Nu/2)*(varTherm + Rx.PD.varShot(Prx, noiseBW)... % thermal + shot
                + Rx.PD.R^2*(Rx.OptAmp.varSigSpont(Prx/Rx.OptAmp.Gain, noiseBW)... % sig-spont
                    + Rx.OptAmp.varSpontSpont(noiseBW, BWopt, Npol))... % spont-spont
                    + varQuantizTx(Prx) + varQuantizRx(Prx)); % quantization noise 
        % Note: Prx is divided by amplifier gain to obtain power at the amplifier input
    else % unamplified system: thermal-noise dominant
        Prx = Tx.Ptx*linkAtt;
        varNoise = 1/ofdm.Nc*ones(1, ofdm.Nu/2)*(varTherm... % thermal noise 
            + Rx.PD.varShot(Prx, noiseBW)... % shot noise 
            + varQuantizTx(Prx) + varQuantizRx(Prx)); % quantization noise 
    end    
    
    % Power allocation
    ofdm.power_allocation(Hch, varNoise, sim.BERtarget, sim.shouldPlot('Power allocation'));
    
%     ber.theory(k) = ofdm.calc_ber(10*log10(ofdm.Pn*Prx/ofdm.dc_bias(ofdm.Pn, Tx.rclip, Tx.Mod.Hdac, Tx.rexdB)./varNoise));
    
    % Montecarlo simulation
    sim.Hch = Hch;
    sim.varNoise = varNoise;
    [ber.count(k), ber.gauss(k), SNRndB, OSNRdB(k)] = ber_dc_ofdm_montecarlo(ofdm, Tx, Fibers, Rx, sim);
    SNRdB(k) = sum(SNRndB.*log2(ofdm.CS))/sum(log2(ofdm.CSn));
    
    if sim.stopSimWhenBERReaches0 && ber.count(k) == 0
        break;
    end 
end

%% Plots
if sim.shouldPlot('BER') && length(ber.count) > 1
    figure(1), hold on, box on
    hline = plot(Tx.PtxdBm, log10(ber.count), '-o', 'LineWidth', 2);
    plot(Tx.PtxdBm, log10(ber.gauss), '-', 'Color', get(hline, 'Color'), 'LineWidth', 2)
%     plot(Tx.PtxdBm, log10(ber.theory), ':k', 'LineWidth', 2)
    legend('Counted', 'Gaussian approximation')
    set(gca, 'FontSize', 12)
    axis([Tx.PtxdBm(1) Tx.PtxdBm(end) -8 0])
    xlabel('Received power (dBm)', 'FontSize', 12)
    ylabel('log_{10}(BER)', 'FontSize', 12)
    set(gca, 'FontSize', 12)
    drawnow
end