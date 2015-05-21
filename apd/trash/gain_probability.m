clear, clc, close all

keff = 0.09;
Gapd = 10;

n = 1:50;

pn = apd_gain_distribution(n, 0.09, 10);

figure
stem(n, pn)
xlabel('n', 'FontSize', 12)
ylabel('p(g_l = n)', 'FontSize', 12)
