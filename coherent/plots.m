%% =============== Figures ===============      
if isfield(sim, 'Plots')
    Plot = @(key) sim.Plots.isKey(key) && sim.Plots(key);
else
    Plot = @(key) false;
end

% Electorde Loss ----------------------------------------------------------
if Plot('Channel frequency response')
    figure(101), hold on
    plot(sim.f/1e9, 20*log10(abs(Tx.Mod.Hel)), '-');
    plot(sim.f/1e9, 20*log10(abs(Tx.filt.H(sim.f/sim.fs))), '-');
    plot(sim.f/1e9, 20*log10(abs(Rx.ADC.filt.H(sim.f/sim.fs))), '-');
    axis([-2*sim.Rs/1e9 2*sim.Rs/1e9 -20 0])
    legend('Modulator', 'Transmitter filter', 'ADC antialiasing')
    xlabel('Frequency (GHz)'); ylabel('Magnitude Response')
    title('Electrode Loss frequency response')
end

% Equalizer frequency response and convergence-----------------------------
if Plot('Equalizer')
    Wx = W{1};
    Wy = W{2};
    AdEq = Rx.AdEq;
    figure(102)
    subplot(221), hold on, box on
    plot(1:length(MSE), MSE(1, :))
    a = axis;
    plot([1 1]*AdEq.Ntrain, [a(3) a(4)], '--k')
    axis([1 length(MSE) a(3:4)])
    xlabel('Iteration'); ylabel('MSE')
    
    subplot(222), hold on, box on
    plot(1:length(MSE), MSE(2, :))
    a = axis;
    plot([1 1]*AdEq.Ntrain, [a(3) a(4)], '--k')
    axis([1 length(MSE) a(3:4)])
    xlabel('Iteration'); ylabel('MSE')
    
    subplot(223), hold on, box on
    for filt = 1:size(Wx, 2)
        [hx, w] = freqz(Wx(:, filt), 1);
        hy = freqz(Wy(:, filt), 1, w);
        hline(1) = plot(w/(2*pi), abs(hx).^2, '-');
        hline(2) = plot(w/(2*pi), abs(hy).^2, '--', 'Color', get(hline(1), 'Color'));
    end
    xlabel('Frequency')
    ylabel('Magnitude')
    legend(hline, {'X pol', 'Y pol'})
    title(sprintf('%s, %d taps', Rx.AdEq.type, Rx.AdEq.Ntaps))
    
    subplot(224), hold on, box on
    for filt = 1:size(Wx, 2)
        [hx, w] = freqz(Wx(:, filt), 1);
        hy = freqz(Wy(:, filt), 1, w);
        hline(1) = plot(w/(2*pi),  unwrap(angle(hx)), '-');
        hline(2) = plot(w/(2*pi),  unwrap(angle(hy)), '--', 'Color', get(hline(1), 'Color'));
    end    
    xlabel('Frequency')
    ylabel('Phase')
    legend('X pol', 'Y pol')
    title(sprintf('%s, %d taps', Rx.AdEq.type, Rx.AdEq.Ntaps))
end

% Differential Group Delay ------------------------------------------------
if Plot('Diff group delay')
    figure(104)
    plot(sim.f/1e9, Fiber.calcDGD(2*pi*sim.f), 'b.');
    axis([-sim.Rs/1e9 sim.Rs/1e9 0 20])
    xlabel('Frequency \omega/2\pi (GHz)'); ylabel('|\Delta\tau(\omega)| (ps)');
    title('Differential group delay vs. \omega');
end

% Eye diagrams ------------------------------------------------------------
if Plot('Eye diagram') 
    eyediagram(abs(Ein(1, sim.Ndiscard:min(sim.Ndiscard+2e3, end))).^2 + ...
    abs(Ein(2, sim.Ndiscard:min(sim.Ndiscard+2e3, end))).^2, 2*sim.Mct)
    title('Eye diagram')
end

% Frequency estimation algorithm ------------------------------------------
 if Plot('Frequency estimation')
     figure(105)
     plot(Df*sim.Rs/1e9)
     xlabel('Iteration')
     ylabel('Frequency (GHz)')
     title('Frequency offset estimation')
 end
 
 % Phase tracker ----------------------------------------------------------
 if Plot('Phase tracker')
     figure(106)
     plot(phiPT.')
     ylabel('Correction phase')
     xlabel('Symbol')
 end
 
 % Phase error plot -------------------------------------------------------
if sim.Plots.isKey('EPLL phase error') && sim.Plots('EPLL phase error')
    figure(204), clf
    subplot(211), hold on, box on
    plot(sim.t, LoopFilterInput)
    xlabel('time (s)')
    ylabel('Loop filter input')       
    subplot(212), hold on, box on
    plot(sim.t, LoopFilterOutput)
    p = polyfit(sim.t, LoopFilterOutput, 1);
    plot(sim.t, polyval(p, sim.t));
    plot(sim.t, sign(p(1))*2*pi*sim.t*Rx.LO.freqOffset)
    legend('VCO phase','Linear fit (for freq offset)',  'Frequency offset ramp')
    xlabel('time (s)')
    ylabel('Phase (rad)')
end
    
    