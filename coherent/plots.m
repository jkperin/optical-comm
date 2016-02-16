%% =============== Figures ===============       
% Electorde Loss ----------------------------------------------------------
if (Plots('ModElectrodeLoss'))
    figure(2)
    plot(sim.f/1e9, 10*log10(abs(Mod.Hel)), '.');
    axis([-Rs/1e9*(M+1)/2 Rs/1e9*(M+1)/2 -5 0])
    xlabel('Frequency (GHz)'); ylabel('Magnitude Response')
    title('Electrode Loss frequency response')
end
% Differential Group Delay ------------------------------------------------
if (Plots('DiffGroupDelay'))
    figure(3)
    plot(sim.f/1e9, Fiber.calcDGD(2*pi*sim.f), 'b.');
    axis([-Rs/1e9*(M+1) Rs/1e9*(M+1) 0 20])
    xlabel('Frequency \omega/2\pi (GHz)'); ylabel('|\Delta\tau(\omega)| (ps)');
    title('Differential group delay vs. \omega');
end