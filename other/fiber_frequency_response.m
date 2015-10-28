%% Fiber small signal frequency response
clear, clc, close all

addpath ../f/


L = [1, 5, 10, 20]*1e3;

fiber = fiber(L(1));

tx.lamb = 1250e-9;
tx.alpha = 2;
f = linspace(0, 200e9, 100);


figure, hold on, grid on, box on
for k= 1:length(L)
    fiber.L = L(k);
    plot(f/1e9, 20*log10(abs(fiber.H(f, tx))))
end

xlabel('Frequency (GHz)')
ylabel('Frequency Response (dB)')
legend('1 km', '5 km', '10 km', '20 km', 'Location', 'SouthWest')
axis([0 200 -10 10])

matlab2tikz('Hfiber.tex')


    
    
