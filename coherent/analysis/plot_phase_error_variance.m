%% plot phase error vs loop filter relaxation frequency for various loop delays
clear, clc, close all

sim.Plots = containers.Map(); 
sim.Plots('Phase error variance') = 0;

csi = sqrt(2)/2;
Delays = [0, 250, 500, 750]*1e-12; 
linewidth = 2*200e3;

sim.Rb = 2*112e9;
sim.M = 4;
sim.Rs = sim.Rb/(2*log2(sim.M));
sim.BERtarget = 1.8e-4;
sim.ModFormat = 'QAM';

figure, hold on, box on
for k = 1:length(Delays)
    [wnOpt, wn, phiError] = optimizePLL(csi, Delays(k), linewidth, sim);
    
    h(k) = plot(wn/1e9, phiError, 'LineWidth', 2);
    plot(wnOpt/1e9, min(phiError), 'o', 'Color', get(h(k), 'Color'), 'MarkerFaceColor', 'w', 'LineWidth', 2)
end
axis([0 1 0 0.03])
xlabel('\omega_n (Grad/s)', 'FontSize', 12)
ylabel('Phase error variance (rad^2)', 'FontSize', 12)
leg = legend(h, {'0 ps', '250 ps (14 symbols)', '500 ps (28 symbols)', '750 ps (42 symbols)'});
set(leg, 'FontSize', 12)
set(gca, 'FontSize', 12)
