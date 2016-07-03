function [ber, ofdm] = dmt_ber(ofdm, Tx, Fibers, Amp, Rx, sim)
%% Calculate BER of pre-amplified IM-DD system
% Inputs:
% - mpam: OFDM class
% - Tx: struct with transmitter paramters
% - Fibers: array of Fiber classes containing the fibers used in
% transmittion i.e., Fibers = [SMF DCF] or Fibers = [SMF];
% - Amp: pre-amplifier using SOA class
% - Rx: struct with receiver parameters
% - sim: struct with simulation parameters

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
    Hmod = @(f) Tx.Mod.H(f); % group delay was already removed
else
    Hmod = @(f) ones(size(f));
end

% Antialiasing filter
Hrx = @(f, fs) Rx.ADC.filt.H(f/fs);

%% Calculate cyclic prefix
Hch = @(f, fs) Hdac(f, fs).*Hmod(f).*Hfiber(f).*Hrx(f, fs);
ofdm.cyclic_prefix(Hch, sim.shouldPlot('Cyclic prefix'));

Hfiber = Hfiber(sim.f);
Hdac = Hdac(sim.f, sim.fs);
Hmod = Hmod(sim.f);
Tx.Mod.H = Hmod;
Hrx = Hrx(sim.f, sim.fs);
Hch = Hch(ofdm.fc, sim.fs); % channel frequency response at the subcarriers frequency

% Noise bandwidth
noiseBW = ofdm.fs/2;
BWopt = Rx.optfilt.fcnorm*sim.fs; % optical filter noise bandwidth; Not divided by 2 because optical filter is a bandpass filter

% Thermal noise
varTherm = Rx.N0*noiseBW; % variance of thermal noise

% Noise std for intensity level Plevel
Npol = 2; % number of polarizations. Npol = 1, if polarizer is present, Npol = 2 otherwise.

% Check if power is at reasonable levels
assert(sqrt(ofdm.var) < 1e-3, 'Power is too high: %.2G (mW)\nCannot improve performance by increasing power.\n', 1e3*sqrt(2*sum(ofdm.Pn)))

%% Calculate BER
Ptx = dBm2Watt(Tx.PtxdBm); % Transmitted power

ber.awgn = zeros(size(Ptx));
ber.count = zeros(size(Ptx));
for k = 1:length(Ptx)
    Tx.Ptx = Ptx(k);
    Tx.Prx = Tx.Ptx*linkAtt;
    
    % Noise variance
    varNoise =  1/ofdm.Nc*ones(1, ofdm.Nu/2)*(varTherm + Rx.PD.varShot(Tx.Prx*Amp.Gain, noiseBW)...
    + Rx.PD.R^2*Amp.varAWGN(Tx.Prx/Amp.Gain, noiseBW, BWopt, Npol));
    % Note: Plevel is divided by amplifier gain to obtain power at the amplifier input

    % Power allocation
    ofdm.power_allocation(Hch, varNoise, sim.BERtarget, sim.shouldPlot('Power allocation'));
    
    % Montecarlo simulation
    sim.Hch = Hch;
    sim.varNoise = varNoise;
    [ber.count(k), ber.awgn(k), ~] = ber_dc_ofdm(ofdm, Tx, Fibers, Amp, Rx, sim)
    
    if sim.stopSimWhenBERReaches0 && ber.count(k) == 0
        break;
    end 
end

%% Plots
if sim.shouldPlot('BER') && length(ber.count) > 1
    figure(1), hold on, box on
    hline = plot(Tx.PtxdBm, log10(ber.count), '-o', 'LineWidth', 2);
    plot(Tx.PtxdBm, log10(ber.awgn), '-', 'Color', get(hline, 'Color'), 'LineWidth', 2)
    legend('Counted', 'Estimated using AWGN approximation')
    set(gca, 'FontSize', 12)
    axis([Tx.PtxdBm(1) Tx.PtxdBm(end) -8 0])
    xlabel('Received power (dBm)', 'FontSize', 12)
    ylabel('BER', 'FontSize', 12)
    set(gca, 'FontSize', 12)
    drawnow
end

if sim.shouldPlot('Channel frequency response')
    fplot = sim.f/1e9;
    figure(107), clf, box on, hold on
    plot(fplot, abs(Hdac).^2)
    plot(fplot, abs(Hmod).^2)
    plot(fplot, abs(Hfiber).^2)
%     plot(fplot, abs(Hpd).^2)
    plot(fplot, abs(Hrx).^2)
    xlabel('Frequency (GHz)', 'FontSize', 12)
    ylabel('|H(f)|^2', 'FontSize', 12)
    leg = legend('DAC', 'Modulator', 'Fiber', 'Receiver filtering');
    set(leg, 'FontSize', 12)
    set(gca, 'FontSize', 12)
    a = axis;
    axis([0 ofdm.fs/1e9 0 a(4)])
    drawnow
end
