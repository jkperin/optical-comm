%% Load experiment results
clear, clc, close all

folder = 'data/waveforms/BER_vs_OSNR_km/'
files = {'PAM4_56G_BER_vs_OSNR_0km_set1',...
            'PAM4_56G_BER_vs_OSNR_0km_set2',...
            'PAM4_56G_BER_vs_OSNR_10km_set1',...
            'PAM4_56G_BER_vs_OSNR_3km_set1',...
            'PAM4_56G_BER_vs_OSNR_3km_set2',...
            'PAM4_56G_BER_vs_OSNR_5km_set1',...
            'PAM4_56G_BER_vs_OSNR_5km_set2'};
    
Experiment = containers.Map();
figure(1), box on, hold on
xlabel('OSNR (dB)', 'FontSize', 12)
ylabel('log_{10}(BER)', 'FontSize', 12)
leg = {};
for k = 1:length(files)
    file = files{k}
    
    s = load([folder file]);
    
    Experiment(file) = struct('OSNRdB', s.experiment.OSNRdB,...
        'BER', s.experiment.BER, 'dacfile', s.experiment.dacfile)
   
    
    figure(1)
    plot(s.experiment.OSNRdB, log10(s.experiment.BER.count), '-o', 'LineWidth', 2)
    leg = [leg file];
end

OSNRdB = linspace(25, 40);
ber_theory = pam_ber_from_osnr(4, OSNRdB, 56e9/4);
plot(OSNRdB, log10(ber_theory), 'k', 'LineWidth', 2)
leg = [leg 'Theory'];
legend(leg)
axis([26 45 -8 0])
save PAM4_experiment_summary Experiment

    