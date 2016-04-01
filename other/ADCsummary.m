%% Plot ADC points
clear, clc, close all

fs = [90
56
46
40
40
40
35
35
35
35
25
25
24
24
24
22
22
20
20
20
16
10.3
10
10
10
10];

ENOB = [5.19
5.68
3.89
3.90
3.03
3.5
3.7
3.2
5
6
4.64
4.00
3.50
2.3
5
3.03
5
4.81
4.60
4.6
4.36
4.96
5.31
5.10
4.74
4.56];

P = [0.667
2
0.381
1.5
3.8
NaN
4.5
4.5
NaN
2
0.088
0.5
1.2
3.8
3.3
3
3
0.0695
10
9
0.435
0.39
0.032
0.24
0.079
0.139];

Tech = [0.032
0.065
0.028
0.065
NaN
NaN
0.18
NaN
NaN
NaN
0.065
0.04
0.09
NaN
0.18
0.13
0.13 
0.032
0.18
0.18
0.065
0.065
0.028
0.04
0.065
0.04];

class = sort(unique(Tech(~isnan(Tech))));

id = (2:length(fs)+1).';
idstr = num2str(id);
idcell = cellstr(idstr);
dx = 0.1; dy1 = 0.1; dy2 = 0.05; % displacement so the text does not overlay the data points
for k = 1:length(class)
    c = class(k);
    idx = ismember(Tech, c);
    f = fs(idx);
    yENOB = ENOB(idx);
    yP = P(idx);
    figure(1), hold on, box on
    h = scatter(f, yENOB, 60, 'filled');
    text(f+dx, yENOB + dy1, idcell(idx))
    hline1(k) = plot(-1, -1, 'o', 'Color', get(h, 'CData'), 'MarkerFaceColor', get(h, 'CData'));
    figure(2), hold on, box on
    scatter(f, yP, 60, 'filled');
    text(f+dx, yP + dy2, idcell(idx))
    hline2(k) = plot(-1, -1, 'o', 'Color', get(h, 'CData'), 'MarkerFaceColor', get(h, 'CData'));
end

figure(1)
xlabel('Sampling Frequency (GHz)', 'FontSize', 12)
ylabel('ENOB (bits)', 'FontSize', 12)
set(gca, 'FontSize', 12)
legend(hline1, cellstr(num2str(class)))
axis([0 90 3 6])

figure(2)
xlabel('Sampling Frequency (GHz)', 'FontSize', 12)
ylabel('Power (W)', 'FontSize', 12)
axis([0 90 0 2])
set(gca, 'FontSize', 12)
legend(hline2, cellstr(num2str(class)))





