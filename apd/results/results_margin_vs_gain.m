%% Compile results from margin_vs_gain.m script
clear, clc, close all

load('margin_vs_gain_results_8-PAM_equally-spaced.mat')

figure, hold on, box on
legs = {};
for n = 1:length(ka)
    hline(n) = plot(BW/1e9, OptMargindB(n, :), '-', 'LineWidth', 2);
    legs = [legs sprintf('ka = %.2f', ka(n))];
end
xlabel('Bandwidth (GHz)')
ylabel('Optimal Margin (dB)')
legend(legs, 'Location', 'SouthEast')
axis([BW([1 end])/1e9 0 9])
grid on
saveas(gca, 'margin_vs_bw', 'emf')

figure, box on, hold on
legs = {};
for n = 1:length(ka)
    plot(BW/1e9, Gopt_margin(n, :), '-', 'Color', get(hline(n), 'Color'), 'LineWidth', 2);
    legs = [legs sprintf('ka = %.2f', ka(n))];
end
xlabel('Bandwidth (GHz)')
ylabel('Optimal Gain')
legend(legs)
axis([BW([1 end])/1e9 0 20])
grid on


% Find gain corresponding to margin pen dB below optimal margin
pendB = 1;

% % MargindB(ka, BW, Gain)
subOptGain = zeros(length(ka), length(BW));
for n = 1:4
    for k = 1:length(BW)
        margin = MargindB(n, k, :);
        margin = margin(:);
        valid = (~isinf(margin) & ~isnan(margin)); % eliminates NaN and Inf
        valid = valid & (Gains.' < Gopt_margin(n, k)); % selects only gains that are below optimal gain
        
        subOptGain(n, k) = interp1(margin(valid), Gains(valid), OptMargindB(n, k) - pendB);
    end
    
    plot(BW/1e9, subOptGain(n, :), '--', 'Color', get(hline(n), 'Color'),  'LineWidth', 2);
end

saveas(gca, 'gain_vs_bw', 'emf')

% figure, hold on, box on
% for n = 1:length(ka)
%     for k = 1:length(BW)
%         margin = MargindB(n, k, :);
%         margin = margin(:);
%         plot(Gains, margin) 
%         drawnow;
%         1;
%     end
% end

