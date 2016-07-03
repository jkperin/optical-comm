function plot_constellation(Y, dataTX, M)
%% Scatterplot with different colors for each symbol
% Inputs:
% - Y: symbols
% - dataTX: symbol identification i.e., integers from 0 to M-1
% - M: constellation order

box on, hold on
for k = 0:M-1
    s = Y(dataTX == k);
    scatter(real(s), imag(s), 10, 'o', 'filled')
end

axis([-1 1 -1 1]*ceil(sqrt(M)));
axis square