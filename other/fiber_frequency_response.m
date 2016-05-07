%% Fiber small signal frequency response
clear, clc, close all

addpath ../f/

m2tikz = matlab2tikz();

L = [1, 5, 10, 80]*1e3;

fiber = fiber(L(1));

tx.lamb = 1550e-9;
tx.alpha = 0;
f = linspace(0, 50e9, 1000);
Colors = {'green', 'blue', 'orange', 'red'};

figure, hold on, grid on, box on
for k= 1:length(L)
    fiber.L = L(k);
    hline(k) = plot(f/1e9, 20*log10(abs(fiber.H(f, tx))));
    plot(f/1e9, 20*log10(abs(fiber.Hlarge_signal(f, tx))), '--', 'Color', get(hline(k), 'Color'))
    % addplot(x, y, line, color, marker, label)
    m2tikz.addplot(f/1e9, 20*log10(abs(fiber.H(f, tx))), '-', Colors{k}, 'none', sprintf('%d km', L(k)/1e3))
end

xlabel('Frequency (GHz)')
ylabel('Frequency Response (dB)')
legend(hline, '1 km', '5 km', '10 km', '20 km', 'Location', 'SouthWest')
axis([0 50 -10 10])

%% Extract curves and formatting from matlab 
% m2tikz.extract(gca)
% m2tikz.write('Hfiber.tex')

%% Extract just formatting from matlab 
m2tikz.extract(gca, 'just axis')
m2tikz.write('Hfiber.tex')

% matlab2tikz('Hfiber.tex')


    
    
