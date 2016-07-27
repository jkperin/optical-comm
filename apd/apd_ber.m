function [ber, mpam, Apd] = apd_ber(mpam, Tx, Fiber, Apd, Rx, sim)
%% Calculate BER of unamplified IM-DD system with APD detector for each transmitted power value in Tx.PtxdBm
% BER is calculated via montecarlo simulation, analytically, AWGN channel, 
% AWGN channel including noise enhancement penalty.

% If equalizer is not defined assume no equalization
if ~isfield(Rx, 'eq')
    Rx.eq.type = 'None';
    Rx.eq.ros = 1; % assuming symbol-rate sampling
end

if strcmpi(Rx.eq.type, 'Fixed TD-SR-LE')
    assert(Rx.eq.ros == 1 && sim.ros.rxDSP == 1, 'If fixed time-domain linear equalizer is used, then oversampling ratio must be 1')
end

% System received pulse shape frequency response
Tx.Ptx = dBm2Watt(Tx.PtxdBm(end));
[~, H] = apd_system_received_pulse_shape(mpam, Tx, Fiber, Apd, Rx, sim); % this is only used if rx.filtering = matched or eq.type = fixed...

% Optimize APD gain
if isfield(sim, 'OptimizeGain') && sim.OptimizeGain
    [Apd.Gain, mpam] = Apd.optGain(mpam, Tx, Fiber, Rx, sim);
    fprintf('Optimal APD Gain = %.2f (%2.f dB)\n', Apd.Gain, Apd.GaindB);
    % if mpam.level_spacing = 'optimized', then Apd.optGain returns mpam 
    % with optimal level spacing
elseif mpam.optimize_level_spacing  %% Level Spacing Optimization
    % Optimize levels using Gaussian approximation
    [~, mpam] = Apd.optimize_PAM_levels(Apd.Gain, mpam, Tx, Fiber, Rx, sim);
    mpam = mpam.norm_levels();
end

%% BER
% Transmitted power
Ptx = dBm2Watt(Tx.PtxdBm);

ber.count = zeros(size(Ptx)); % counted BER
ber.enum = zeros(size(Ptx)); % analysis assuming Gaussian stats
ber.awgn = zeros(size(Ptx)); % AWGN approximation (includes noise enhancement penalty)
run_montecarlo = true;
for k = 1:length(Ptx)
    Tx.Ptx = Ptx(k);
            
    % Montecarlo simulation
    if run_montecarlo
        [ber.count(k), mpamOpt] = ber_apd_montecarlo(mpam, Tx, Fiber, Apd, Rx, sim);
        sim.mpamOpt = mpamOpt; % optimized by sweeping thresholds if levels are equally spaced
    end
    
    % BER using Gaussian stats approximation for shot noise (enumeration)
    ber.enum(k) = ber_apd_enumeration(mpam, Tx, Fiber, Apd, Rx, sim);
    
    % BER using AWGN system approximation including noise enhacement
    ber.awgn(k) = ber_apd_awgn(mpam, Tx, Fiber, Apd, Rx, sim);
    
    if isfield(sim, 'terminateWhenBERReaches0') && sim.terminateWhenBERReaches0 && ber.count(k) == 0
        run_montecarlo = false;
        sim.mpamOpt = []; % discard optimized thresholds. Now ber_apd_enumeration will use thresholds at the mid point
    end
end

%% Plot BER
if sim.shouldPlot('BER') && length(ber.count) > 1
    PrxdBm = Tx.PtxdBm - 10*log10(Fiber.link_attenuation(Tx.Laser.wavelength));
    figure(1), hold on, box on
    hline = plot(PrxdBm, log10(ber.count), 'o');
    plot(PrxdBm, log10(ber.enum), '-', 'Color', get(hline, 'Color'))
    plot(PrxdBm, log10(ber.awgn), '--', 'Color', get(hline, 'Color'))
    legend('Counted', 'Enumeration', 'AWGN approximation',...
        'Location', 'SouthWest')
    xlabel('Received Power (dBm)')
    ylabel('log_{10}(BER)')
    grid on
    axis([Tx.PtxdBm(1) Tx.PtxdBm(end) -8 0])  
    drawnow
end

% Frequency respnose
if sim.shouldPlot('Frequency Response')   
    figure(100), clf, hold on, box on 
    plot(sim.f/1e9, abs(H.dac).^2)
    plot(sim.f/1e9, abs(H.mod).^2)
    plot(sim.f/1e9, abs(H.fiber).^2)
    leg = {'DAC', 'Modulator', 'Fiber'};
        
    if ~isinf(Apd.BW)
        plot(sim.f/1e9, abs(H.apd).^2)
        leg = [leg 'APD'];
    end
    
    if ~isinf(Apd.BW) && sim.WhiteningFilter
        plot(sim.f/1e9, abs(H.w).^2)
        leg = [leg 'Noise Whitening'];
    end
    
    plot(sim.f/1e9, abs(H.rx).^2)
    leg = [leg 'Receiver filter'];       
    xlabel('Frequency (GHz')
    ylabel('Frequency Response')
    a = axis;
    axis([0 2*mpam.Rs/1e9 a(3) a(4)])
    legend(leg)
end