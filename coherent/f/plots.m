%% =============== Figures ===============      
if isfield(sim, 'Plots')
    Plot = @(key) sim.Plots.isKey(key) && sim.Plots(key);
else
    Plot = @(key) false;
end

% Electorde Loss ----------------------------------------------------------
if Plot('Channel frequency response')
    figure(101), hold on
    plot(sim.f/1e9, 20*log10(abs(Tx.Mod.H)), '-');
    plot(sim.f/1e9, 20*log10(abs(Tx.filt.H(sim.f/sim.fs))), '-');
    plot(sim.f/1e9, 20*log10(abs(Rx.ADC.filt.H(sim.f/sim.fs))), '-');
    axis([-2*sim.Rs/1e9 2*sim.Rs/1e9 -20 0])
    legend('Modulator', 'Transmitter filter', 'ADC antialiasing')
    xlabel('Frequency (GHz)'); ylabel('Magnitude Response')
    title('Electrode Loss frequency response')
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
    title('Eye diagram: transmitted electric field')
    
    eyediagram(abs(Erec(1, sim.Ndiscard:min(sim.Ndiscard+2e3, end))).^2 + ...
    abs(Erec(2, sim.Ndiscard:min(sim.Ndiscard+2e3, end))).^2, 2*sim.Mct)
    title('Eye diagram: received electric field')
end

% Frequency estimation algorithm ------------------------------------------
 if Plot('Frequency estimation')
     figure(105)
     plot(Df*sim.Rs/1e9)
     xlabel('Iteration')
     ylabel('Frequency (GHz)')
     title('Frequency offset estimation')
 end
 

 

    
    