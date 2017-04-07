function plot_ber_struct(x, ber)
names = fieldnames(ber);
hold on, box on
idx = [];
for k = 1:length(names)   
    hline(k) = plot(x, log10(getfield(ber, names{k})), '-o', 'LineWidth', 2, 'DisplayName', names{k});
end
axis([x(1) x(end) -8 0])
xlabel('Received power (dBm)', 'FontSize', 12)
ylabel('log_{10}(BER)', 'FontSize', 12)
set(gca, 'FontSize', 12)
legend('-dynamiclegend')
drawnow