%% Process levels plots
clear, clc

addpath ../../f/

figure(1)
h = allchild(gca);

n = 1;
m = 1;
p = 1;
c = [];
for k = 1:length(h)
    if strcmpi(class(h(k)), 'matlab.graphics.chart.primitive.Line')
        if length(h(k).XData) == 2
            t(n) = h(k).XData(1);
            n = n + 1;
            continue
        end
        
        c(m).x = h(k).XData;
        c(m).y = h(k).YData;
%         f = fit(c(m).x.',c(m).y.','gauss1');
        c(m).xnew = linspace(c(m).x(1), c(m).x(end), 1e3); 
        c(m).ynew = interp1(c(m).x, c(m).y, c(m).xnew);
%         c(m).fitted = f(c(m).xnew);
        m = m + 1;
    elseif strcmpi(class(h(k)), 'matlab.graphics.primitive.Text')
        l(p) = h(k).Position(1);
        p = p + 1;
    end
end

xmax = max([c(1).x c(2).x c(3).x c(4).x])
ymax = max([c(1).ynew c(2).ynew c(3).ynew c(4).ynew])

xmax = 8e-4;
ymax = 1.5e5;

% Normalize everything and plot
levels = l/xmax;
thresholds = t/xmax

figure, hold on, box on
for k = 1:4
    plot(c(k).xnew/xmax, c(k).ynew/ymax, 'k', 'LineWidth', 2)
    axis tight
end

m = matlab2tikz(gca);
m.write('SOA-equally_spaced.tex')



