%% Chirped fiber IM-DD frequency response
clear, clc, close all

f = linspace(0, 30, 1e3)*1e9;
lambda = 1550e-9;
alpha = 2;

Fiber = fiber(10e3);

figure, hold on, box on
plot(f/1e9, abs(Fiber.Himdd(f, lambda, 0, 'large signal')).^2)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -1, 'large signal')).^2)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -2, 'large signal')).^2)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -3, 'large signal')).^2)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -4, 'large signal')).^2)
xlabel('Frequency (GHz)', 'FontSize', 12)
ylabel('|H_{IM-DD}(f)|^2', 'FontSize', 12)
leg = legend('\alpha = 0', '-1', '-2', '-3', '-4');
set(leg, 'FontSize', 12);
set(gca, 'FontSize', 12)
title(sprintf('Large-signal IM-DD frequency response for L = %.2f km', Fiber.L/1e3))

figure, hold on, box on
plot(f/1e9, abs(Fiber.Himdd(f, lambda, 0, 'small signal')).^2)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -1, 'small signal')).^2)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -2, 'small signal')).^2)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -3, 'small signal')).^2)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -4, 'small signal')).^2)
xlabel('Frequency (GHz)', 'FontSize', 12)
ylabel('|H_{IM-DD}(f)|^2', 'FontSize', 12)
leg = legend('\alpha = 0', '-1', '-2', '-3', '-4');
set(leg, 'FontSize', 12);
set(gca, 'FontSize', 12)
title(sprintf('Small-signal IM-DD frequency response for L = %.2f km', Fiber.L/1e3))