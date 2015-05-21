clear, clc, close all
% Plot curves

% 1. ENOB = 5, preemphasis, no shot no RIN
% 2. ENOB = 6, preemphais, no shot no RIN
% 3. ENOB = 5, palloc, no shot no RIN
% 4. ENOB = 6, palloc, no shot no RIN
% 5. ENOB = 5, preemphasis, shot & RIN
% 6. ENOB = 6, preemphais, shot & RIN
% 7. ENOB = 5, palloc, shot & RIN
% 8. ENOB = 6, palloc, shot & RIN

% preemphasis = circle
% palloc = square
% no shot no rin = dashed
% rin & shot = solid

load('data/results_CS_16QAM')

L = 0:1:5;

figure, hold on, box on, grid on
% ENOB = 5 
plot(-10, -10, '-ok', 'LineWidth', 1.5)
plot(-10, -10, '-sk', 'LineWidth', 1.5)
plot(-10, -10, '--ok', 'LineWidth', 1.5)
plot(-10, -10, '--sk', 'LineWidth', 1.5)
plot(L, results_CS_16QAM{1}.power_pen_ook_m, '--o', 'Color', 'b', 'LineWidth', 1.5)
plot(L, results_CS_16QAM{3}.power_pen_ook_m, '--s', 'Color', 'b', 'LineWidth', 1.5)
plot(L, results_CS_16QAM{5}.power_pen_ook_m, '-o', 'Color', 'b', 'LineWidth', 1.5)
plot(L, results_CS_16QAM{7}.power_pen_ook_m, '-s', 'Color', 'b', 'LineWidth', 1.5)

% % ENOB = 6
% plot(L, results_CS_16QAM{2}.power_pen_ook_m, '--o', 'Color', 'r', 'LineWidth', 1.5)
% plot(L, results_CS_16QAM{4}.power_pen_ook_m, '--s', 'Color', 'r', 'LineWidth', 1.5)
% plot(L, results_CS_16QAM{6}.power_pen_ook_m, '-o', 'Color', 'r', 'LineWidth', 1.5)
% plot(L, results_CS_16QAM{8}.power_pen_ook_m, '-s', 'Color', 'r', 'LineWidth', 1.5)

xlabel('Fiber Length (km)', 'FontSize', 12)
ylabel('Optical Power Penalty (dB)', 'FontSize', 12)
axis([0 5 0 19])
set(gca, 'xtick', L)
legend('Including Shot & RIN | Const. bit loading & preemphasis', 'Including Shot & RIN | Opt. bit loading & power allocation',...
    'Const. bit loading & preemphasis', 'Opt. bit loading & power allocation', 'Location', 'SouthEast');



% 1. ENOB = 6, preemphais, no shot no RIN
% 2. ENOB = 6, palloc, no shot no RIN
% 3. ENOB = 6, preemphais, shot & RIN
% 4. ENOB = 6, palloc, shot & RIN

% preemphasis = circle
% palloc = square
% no shot no rin = dashed
% rin & shot = solid

load('data/results_CS_64QAM')

L = 0:1:5;

% ENOB = 6
plot(L, results_CS_64QAM{1}.power_pen_ook_m, '--o', 'Color', 'r', 'LineWidth', 1.5)
plot(L, results_CS_64QAM{2}.power_pen_ook_m, '--s', 'Color', 'r', 'LineWidth', 1.5)
plot(L, results_CS_64QAM{3}.power_pen_ook_m, '-o', 'Color', 'r', 'LineWidth', 1.5)
plot(L, results_CS_64QAM{4}.power_pen_ook_m, '-s', 'Color', 'r', 'LineWidth', 1.5)

xlabel('Fiber Length (km)', 'FontSize', 12)
ylabel('Optical Power Penalty (dB)', 'FontSize', 12)
axis([0 5 0 20])
set(gca, 'xtick', L)

matlab2tikz('tikz/pp_vs_L.tikz')
