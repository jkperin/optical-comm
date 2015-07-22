clear, clc, close all

load('soa_awgn_results_req_Ptx_vs_Fn.mat')

pin.Ptxreq(1).eq_spaced = -13.696431025987138;
pin.Ptxreq(1).optimized = -13.713654564942999;
pin.Ptxreq(2).eq_spaced = -10.851632818179333;
pin.Ptxreq(2).optimized = -10.899533093257471;
pin.Ptxreq(3).eq_spaced = -7.947794653274341;
pin.Ptxreq(3).optimized = -8.131931360279786;

%% Power improvement
figure, grid on, hold on, box on
Colors = {'b', 'r', 'g', 'k'};
plot(Fn, pin.Ptxreq(1).eq_spaced - PtxdBm_req{1}.eq_spaced, 'LineWidth', 1.5, 'Marker', 'o', 'Color', Colors{1});
plot(Fn, pin.Ptxreq(2).eq_spaced - PtxdBm_req{2}.eq_spaced, 'LineWidth', 1.5, 'Marker', 'o', 'Color', Colors{2});
plot(Fn, pin.Ptxreq(3).eq_spaced - PtxdBm_req{3}.eq_spaced, 'LineWidth', 1.5, 'Marker', 'o', 'Color', Colors{3});

plot(Fn, pin.Ptxreq(1).eq_spaced - PtxdBm_req{1}.optimized, 'LineWidth', 1.5, 'Marker', 's', 'Color', Colors{1});
plot(Fn, pin.Ptxreq(2).eq_spaced - PtxdBm_req{2}.optimized, 'LineWidth', 1.5, 'Marker', 's', 'Color', Colors{2});
plot(Fn, pin.Ptxreq(3).eq_spaced - PtxdBm_req{3}.optimized, 'LineWidth', 1.5, 'Marker', 's', 'Color', Colors{3});

xlabel('Noise Figure (dB)', 'Fontsize', 12)
ylabel('Optical Power Margin Improvement (dB)', 'Fontsize', 12)
legend('4-PAM', '8-PAM', '16-PAM')

figure, grid on, hold on, box on
Colors = {'b', 'r', 'g', 'k'};
plot(Fn, pin.Ptxreq(1).eq_spaced - PtxdBm_req_polarizer{1}.eq_spaced, 'LineWidth', 1.5, 'Marker', 'o', 'Color', Colors{1});
plot(Fn, pin.Ptxreq(2).eq_spaced - PtxdBm_req_polarizer{2}.eq_spaced, 'LineWidth', 1.5, 'Marker', 'o', 'Color', Colors{2});
plot(Fn, pin.Ptxreq(3).eq_spaced - PtxdBm_req_polarizer{3}.eq_spaced, 'LineWidth', 1.5, 'Marker', 'o', 'Color', Colors{3});

plot(Fn, pin.Ptxreq(1).eq_spaced - PtxdBm_req_polarizer{1}.optimized, 'LineWidth', 1.5, 'Marker', 's', 'Color', Colors{1});
plot(Fn, pin.Ptxreq(2).eq_spaced - PtxdBm_req_polarizer{2}.optimized, 'LineWidth', 1.5, 'Marker', 's', 'Color', Colors{2});
plot(Fn, pin.Ptxreq(3).eq_spaced - PtxdBm_req_polarizer{3}.optimized, 'LineWidth', 1.5, 'Marker', 's', 'Color', Colors{3});

xlabel('Noise Figure (dB)', 'Fontsize', 12)
ylabel('Optical Power Margin Improvement (dB)', 'Fontsize', 12)
legend('4-PAM', '8-PAM', '16-PAM')