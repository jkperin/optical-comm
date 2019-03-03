
S264 = load('results/capacity_vs_pump_power_EDF=corning high NA_pump=60mW_980nm_ChDf=50GHz_L=264_x_50km.mat');
S220 = load('results/edf_6m_1dB_margin/capacity_vs_pump_power_EDF=corning high NA_pump=60mW_980nm_ChDf=50GHz_L=220_x_50km.mat');
S = load('recirculating_loop_simulation_Pump=60mW_fix_ideal_GFF.mat')

lnm = S220.Signal.lnm;

figure(800), hold on, box on
plot(lnm, S220.nlin.S.PdBm, 'LineWidth', 2)
plot(lnm, S264.nlin.S.PdBm, 'LineWidth', 2)
xlabel('Wavelength (nm)', 'FontSize', 12)
ylabel('Optimal power profile (dBm)', 'FontSize', 12)
legend('220 spans', '264 spans', 'Location', 'SouthEast')
ylim([-20 -15])
set(gca, 'FontSize', 12)

figure(801), hold on, box on
plot(lnm, S220.lin.num.SE, 'LineWidth', 2)
plot(lnm, S220.nlin.num.SE, 'LineWidth', 2)
plot(lnm, S264.nlin.num.SE, 'LineWidth', 2)
xlabel('Wavelength (nm)', 'FontSize', 12)
ylabel('Spectral efficiency (bit/s/Hz)', 'FontSize', 12)
legend('220 spans', '264 spans', 'Location', 'SouthEast')
set(gca, 'FontSize', 12)
ylim([1.5 5.5])

figure(802), hold on, box on
plot(lnm, FixGFFperAmp, 'LineWidth', 2)
xlabel('Wavelength (nm)', 'FontSize', 12)
ylabel('GFF profile (dB)', 'FontSize', 12)
set(gca, 'FontSize', 12)

figure(803), hold on, box on
plot(lnm, S.SE, 'LineWidth', 2)
plot(lnm, S264.nlin.num.SE, 'LineWidth', 2)
xlabel('Wavelength (nm)', 'FontSize', 12)
ylabel('Spectral efficiency (bit/s/Hz)', 'FontSize', 12)
set(gca, 'FontSize', 12)
legend('Recirculating loop simulation (220 spans)', 'Optimization (264 spans)', 'Location', 'SouthEast')

