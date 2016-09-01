function [ber, OSNRdB, mpam] = preamplified_sys_ber(mpam, Tx, Fibers, Amp, Rx, sim)
%% Calculate BER of pre-amplified IM-DD system
% Inputs:
% - mpam: PAM class
% - Tx: struct with transmitter paramters
% - Fibers: array of Fiber classes containing the fibers used in
% transmittion i.e., Fibers = [SMF DCF] or Fibers = [SMF];
% - Amp: pre-amplifier using SOA class
% - Rx: struct with receiver parameters
% - sim: struct with simulation parameters
 
%% Calculate total (residual) dispersion
Dtotal = 0;
for k = 1:length(Fibers)
    fiberk = Fibers(k);
    
    Dtotal = Dtotal + fiberk.D(Tx.Laser.lambda)*fiberk.L;
end

fprintf('Total dispersion at %.2f nm: %.3f ps/nm\n', Tx.Laser.lambda*1e9, 1e3*Dtotal)

% Create fiber object with equivalent dispersion.
equivFiber = fiber(1, @(lamb) 0, @(lamb) Dtotal);
Hfiber = equivFiber.Himdd(sim.f, Tx.Laser.wavelength, Tx.Mod.alpha, 'large signal'); % frequency response of the channel (used in designing the equalizer)

%% Calculate components frequency response
f = sim.f/sim.fs;

% pulse shape
Hpshape = mpam.Hpshape(sim.f);

% ZOH and DAC frequency response
Nhold = sim.Mct/Tx.DAC.ros;
hZOH = 1/Nhold*ones(1, Nhold);
Hdac = Tx.DAC.filt.H(f).*freqz(hZOH, 1, f);

% Modulator frequency response
if isfield(Tx, 'Mod')
    Hmod = Tx.Mod.H; % group delay was already removed
else
    Hmod = ones(size(f));
end

% Antialiasing filter
Hrx = Rx.ADC.filt.H(f);

%% Calculate BER
Ptx = dBm2Watt(Tx.PtxdBm); % Transmitted power

ber.gauss = zeros(size(Ptx));
ber.count = zeros(size(Ptx));
OSNRdB = zeros(size(Ptx));
for k = 1:length(Ptx)
    Tx.Ptx = Ptx(k);
         
    % Montecarlo simulation
    [ber.count(k), ber.gauss(k), OSNRdBstruct] = ber_preamp_sys_montecarlo(mpam, Tx, Fibers, Amp, Rx, sim);
    OSNRdB(k) = OSNRdBstruct.theory;
    
    if sim.stopSimWhenBERReaches0 && ber.count(k) == 0
        break;
    end 
end

%% Plots
if sim.shouldPlot('BER') && length(ber.count) > 1
    figure(1), hold on, box on
    hline(1) = plot(OSNRdB, log10(ber.count), '-o', 'LineWidth', 2, 'DisplayName', 'Counted');
    hline(2) = plot(OSNRdB, log10(ber.gauss), '-', 'Color', get(hline(1), 'Color'), 'LineWidth', 2, 'DisplayName', 'Gaussian Approximation');
    hline(3) = plot(OSNRdB(OSNRdB ~= 0), log10(pam_ber_from_osnr(mpam.M, OSNRdB(OSNRdB ~= 0), mpam.Rs/2, -Inf)), ':k', 'LineWidth', 2, 'DisplayName', 'Sig-spont limit');
    legend('-DynamicLegend')
    set(gca, 'FontSize', 12)
    axis([OSNRdB(1) OSNRdB(find(OSNRdB ~= 0, 1, 'last')) -8 0])
    xlabel('OSNR (dB)', 'FontSize', 12)
    ylabel('log_{10}(BER)', 'FontSize', 12)
    set(gca, 'FontSize', 12)
    drawnow        
%     figure(1), hold on, box on
%     hline(1) = plot(Tx.PtxdBm, log10(ber.count), '-o', 'LineWidth', 2);
%     hline(2) = plot(Tx.PtxdBm, log10(ber.awgn), '-', 'Color', get(hline(1), 'Color'), 'LineWidth', 2);
%     legend(hline, {'Counted', 'AWGN approximation'})
%     set(gca, 'FontSize', 12)
%     axis([Tx.PtxdBm(1) Tx.PtxdBm(end) -8 0])
%     xlabel('Received power (dBm)', 'FontSize', 12)
%     ylabel('log_{10}(BER)', 'FontSize', 12)
%     set(gca, 'FontSize', 12)
%     drawnow
end

if sim.shouldPlot('Channel frequency response')
    fplot = sim.f/1e9;
    figure(107), clf, box on, hold on
    plot(fplot, abs(Hpshape).^2)
    plot(fplot, abs(Hdac).^2)
    plot(fplot, abs(Hmod).^2)
    plot(fplot, abs(Hfiber).^2)
%     plot(fplot, abs(Hpd).^2)
    plot(fplot, abs(Hrx).^2)
    xlabel('Frequency (GHz)', 'FontSize', 12)
    ylabel('|H(f)|^2', 'FontSize', 12)
    leg = legend('Pulse Shape', 'DAC', 'Modulator', 'Fiber', 'Receiver filtering');
    set(leg, 'FontSize', 12)
    set(gca, 'FontSize', 12)
    a = axis;
    axis([0 2*mpam.Rs/1e9 0 a(4)])
    drawnow
end
    