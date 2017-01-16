%% Evaluate effect of flicker noise on phase error variance
clear, clc, close all

addpath ../f/
addpath ../../f/

Ncpr = 1;
csi = sqrt(2)/2;
Delays = [0, 250, 500, 750]*1e-12; 
totalLinewidth = 2*500e3;

% ka = s1/2 of equation (3) of 
% "Differential carrier phase recovery for QPSK optical coherent systems with integrated tunable lasers"
% DOI: 10.1364/OE.21.010166
% ka = 0.75e10; % DFB
ka = 1.7e12; % Integrated tunable laser: digital supermode-distributed Bragg reflector (DS-DBR)

sim.Rb = 2*112e9;
sim.ModFormat = QAM(4, sim.Rb/2);
sim.Rs = sim.ModFormat.Rs; %sim.Rb/(2*log2(sim.M));
sim.BERtarget = 1e-4;
SNRdB = 11;

% wn = linspace(0, 2*pi*1, 100)*1e9;
wn = 2*pi*1e9*(0:0.01:0.3);
for k = 1:length(Delays)
    [wnOpt, minVarPhiError] = optimizePLL(csi, Delays(k), totalLinewidth, Ncpr, sim);
    
    [varPhiError, nPN, nFlicker, nAWGN] = phase_error_variance_with_flicker(ka, csi, wn, Ncpr, Delays(k), totalLinewidth, SNRdB, sim.Rs, true);
    
    figure(2), hold on
    plot(wn, nPN, wn, nAWGN)
    
    figure(1), hold on, box on
    h(k) = plot(wn/1e9, varPhiError, 'k', 'LineWidth', 2);
    plot(wnOpt/1e9, minVarPhiError, 'ko', 'Color', get(h(k), 'Color'), 'MarkerFaceColor', 'w', 'LineWidth', 2)
    drawnow
end
axis([0 1 0 0.03])
xlabel('\omega_n (Grad/s)', 'FontSize', 12)
ylabel('Phase error variance (rad^2)', 'FontSize', 12)
leg = legend(h, {'0 ps', '250 ps (14 symbols)', '500 ps (28 symbols)', '750 ps (42 symbols)'});
set(leg, 'FontSize', 12)
set(gca, 'FontSize', 12)
