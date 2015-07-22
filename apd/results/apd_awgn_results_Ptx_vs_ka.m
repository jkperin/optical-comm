clear, clc, close all

load('apd_awgn_results_Ptx_vs_ka.mat')

pin.Ptxreq(1).eq_spaced = -13.696431025987138;
pin.Ptxreq(1).optimized = -13.713654564942999;
pin.Ptxreq(2).eq_spaced = -10.851632818179333;
pin.Ptxreq(2).optimized = -10.899533093257471;
pin.Ptxreq(3).eq_spaced = -7.947794653274341;
pin.Ptxreq(3).optimized = -8.131931360279786;

%% Optimal Gain
figure, grid on, hold on, box on
hplot(1) = plot(ka, Gapd_opt{1}.eq_spaced, 'LineWidth', 1.5);
hplot(2) = plot(ka, Gapd_opt{2}.eq_spaced, 'LineWidth', 1.5);
hplot(3) = plot(ka, Gapd_opt{3}.eq_spaced, 'LineWidth', 1.5);

plot(ka, Gapd_opt{1}.optimized, 'LineWidth', 1.5, 'LineStyle', '--', 'Color', get(hplot(1), 'Color'))
plot(ka, Gapd_opt{2}.optimized, 'LineWidth', 1.5, 'LineStyle', '--', 'Color', get(hplot(2), 'Color'))
plot(ka, Gapd_opt{3}.optimized, 'LineWidth', 1.5, 'LineStyle', '--', 'Color', get(hplot(3), 'Color'))
xlabel('k_a', 'Fontsize', 12)
ylabel('Optimal Gain (linear units)', 'Fontsize', 12)
legend('4-PAM', '8-PAM', '16-PAM')
axis([0.05 1 0 25])

%% Power improvement
figure, grid on, hold on, box on
Colors = {'b', 'r', 'g', 'k'};
m = 1;
for k = [1 6 11 12]
    plot(ka, pin.Ptxreq(1).eq_spaced - PtxdBm_req{1}(k).eq_spaced, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', Colors{m});
    m = m + 1;
end
m = 1;
for k = [1 6 11 12]
    plot(ka, pin.Ptxreq(1).eq_spaced - PtxdBm_req{1}(k).optimized, 'LineWidth', 1.5, 'LineStyle', '--', 'Color', Colors{m});
    m = m + 1;
end
xlabel('k_a', 'Fontsize', 12)
ylabel('Optical Power Margin Improvement (dB)', 'Fontsize', 12)
title('4-PAM')
legend('G = 5', 'G = 10', 'G = 15', 'Optimal')

figure, grid on, hold on, box on
Colors = {'b', 'r', 'g', 'k'};
m = 1;
for k = [1 6 11 12]
    plot(ka, pin.Ptxreq(2).eq_spaced - PtxdBm_req{2}(k).eq_spaced, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', Colors{m});
    m = m + 1;
end
m = 1;
for k = [1 6 11 12]
    plot(ka, pin.Ptxreq(2).eq_spaced - PtxdBm_req{2}(k).optimized, 'LineWidth', 1.5, 'LineStyle', '--', 'Color', Colors{m});
    m = m + 1;
end
xlabel('k_a', 'Fontsize', 12)
ylabel('Optical Power Margin Improvement (dB)', 'Fontsize', 12)
title('8-PAM')
legend('G = 5', 'G = 10', 'G = 15', 'Optimal')

figure, grid on, hold on, box on
Colors = {'b', 'r', 'g', 'k'};
m = 1;
for k = [1 6 11 12]
    plot(ka, pin.Ptxreq(3).eq_spaced - PtxdBm_req{3}(k).eq_spaced, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', Colors{m});
    m = m + 1;
end
m = 1;
for k = [1 6 11 12]
    plot(ka, pin.Ptxreq(3).eq_spaced - PtxdBm_req{3}(k).optimized, 'LineWidth', 1.5, 'LineStyle', '--', 'Color', Colors{m});
    m = m + 1;
end
xlabel('k_a', 'Fontsize', 12)
ylabel('Optical Power Margin Improvement (dB)', 'Fontsize', 12)
title('16-PAM')
legend('G = 5', 'G = 10', 'G = 15', 'Optimal')

% %% Required Power
% figure, grid on, hold on, box on
% leg = {};
% for k = 1:length(Gapd)
%     hplot = plot(ka, PtxdBm_req{1}(k).eq_spaced, 'LineWidth', 1.5);
%     plot(ka, PtxdBm_req{1}(k).optimized, 'LineWidth', 1.5, 'LineStyle', '--', 'Color', get(hplot, 'Color'));
% %     leg = [leg sprintf('G = %d (eq)', Gapd(k)) sprintf('G = %d (opt)', Gapd(k))];
% end
% xlabel('k_a', 'Fontsize', 12)
% ylabel('Required Transmitted Optical Power (dBm)', 'Fontsize', 12)
% legend(leg)
% title('4-PAM')
% 
% figure, grid on, hold on, box on
% leg = {};
% for k = 1:length(Gapd)
%     hplot = plot(ka, PtxdBm_req{2}(k).eq_spaced, 'LineWidth', 1.5);
%     plot(ka, PtxdBm_req{2}(k).optimized, 'LineWidth', 1.5, 'LineStyle', '--', 'Color', get(hplot, 'Color'));
%     leg = [leg sprintf('G = %d (eq)', Gapd(k)) sprintf('G = %d (opt)', Gapd(k))];
% end
% xlabel('k_a', 'Fontsize', 12)
% ylabel('Required Transmitted Optical Power (dBm)', 'Fontsize', 12)
% legend(leg)
% title('8-PAM')
% 
% figure, grid on, hold on, box on
% leg = {};
% for k = 1:length(Gapd)
%     hplot = plot(ka, PtxdBm_req{3}(k).eq_spaced, 'LineWidth', 1.5);
%     plot(ka, PtxdBm_req{3}(k).optimized, 'LineWidth', 1.5, 'LineStyle', '--', 'Color', get(hplot, 'Color'));
%     leg = [leg sprintf('G = %d (eq)', Gapd(k)) sprintf('G = %d (opt)', Gapd(k))];
% end
% xlabel('k_a', 'Fontsize', 12)
% ylabel('Required Transmitted Optical Power (dBm)', 'Fontsize', 12)
% legend(leg)
% title('16-PAM')
