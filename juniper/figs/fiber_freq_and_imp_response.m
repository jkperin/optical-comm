%% Chirped fiber IM-DD frequency response
clear, clc, close all

addpath ../../f

f = linspace(0, 50, 1e3)*1e9;
lambda = 1550e-9;
alpha = 2;

Fiber = fiber(10e3);
%% Large signal, fixed length, variable chirp
figure, hold on, box on
plot(f/1e9, abs(Fiber.Himdd(f, lambda, 0, 'large signal')).^2, 'LineWidth', 4)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -1, 'large signal')).^2, 'LineWidth', 4)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -2, 'large signal')).^2, 'LineWidth', 4)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -3, 'large signal')).^2, 'LineWidth', 4)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -4, 'large signal')).^2, 'LineWidth', 4)
xlabel('Frequency (GHz)', 'FontSize', 18)
ylabel('|H_{IM-DD}(f)|^2', 'FontSize', 18)
leg = legend('\alpha = 0', '-1', '-2', '-3', '-4');
set(leg, 'FontSize', 18);
set(gca, 'FontSize', 18)
% title(sprintf('Large-signal IM-DD frequency response for L = %.2f km', Fiber.L/1e3))
saveas(gca, sprintf('Hfiber_large_signal_L=%dkm.png', Fiber.L/1e3))

%% Small signal, fixed length, variable chirp
figure, hold on, box on
plot(f/1e9, abs(Fiber.Himdd(f, lambda, 0, 'small signal')).^2, 'LineWidth', 4)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -1, 'small signal')).^2, 'LineWidth', 4)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -2, 'small signal')).^2, 'LineWidth', 4)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -3, 'small signal')).^2, 'LineWidth', 4)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -4, 'small signal')).^2, 'LineWidth', 4)
xlabel('Frequency (GHz)', 'FontSize', 18)
ylabel('|H_{IM-DD}(f)|^2', 'FontSize', 18)
leg = legend('\alpha = 0', '-1', '-2', '-3', '-4');
set(leg, 'FontSize', 18);
set(gca, 'FontSize', 18)
% title(sprintf('Small-signal IM-DD frequency response for L = %.2f km', Fiber.L/1e3))
saveas(gca, sprintf('Hfiber_small_signal_L=%dkm.png', Fiber.L/1e3))

%% Large signal, variable length, fixed chirp
figure, hold on, box on
alpha = 0;
Fiber.L = 1e3;
plot(f/1e9, abs(Fiber.Himdd(f, lambda, alpha, 'large signal')).^2, 'LineWidth', 4, 'DisplayName', sprintf('L = %d km', Fiber.L/1e3))
Fiber.L = 5e3;
plot(f/1e9, abs(Fiber.Himdd(f, lambda, alpha, 'large signal')).^2, 'LineWidth', 4, 'DisplayName', sprintf('L = %d km', Fiber.L/1e3))
Fiber.L = 10e3;
plot(f/1e9, abs(Fiber.Himdd(f, lambda, alpha, 'large signal')).^2, 'LineWidth', 4, 'DisplayName', sprintf('L = %d km', Fiber.L/1e3))
Fiber.L = 15e3;
plot(f/1e9, abs(Fiber.Himdd(f, lambda, alpha, 'large signal')).^2, 'LineWidth', 4, 'DisplayName', sprintf('L = %d km', Fiber.L/1e3))
Fiber.L = 20e3;
plot(f/1e9, abs(Fiber.Himdd(f, lambda, alpha, 'large signal')).^2, 'LineWidth', 4, 'DisplayName', sprintf('L = %d km', Fiber.L/1e3))
xlabel('Frequency (GHz)', 'FontSize', 18)
ylabel('|H_{IM-DD}(f)|^2', 'FontSize', 18)
legend('-DynamicLegend')
set(leg, 'FontSize', 18);
set(gca, 'FontSize', 18)
% title(sprintf('Large-signal IM-DD frequency response for alpha = %.1f km', alpha))
saveas(gca, sprintf('Hfiber_large_signal_alpha=%dkm.png', alpha))

%% First notch location
Lkm = 0:30;
Alpha = 0:-1:-3;
fpos = f(f >= 0);
figure, box on, hold on
% fnotch(k) = 1/(2*pi)*sqrt(-pi/(beta2*Fiber.L)); % for alpha = 0
for kk = 1:length(Alpha)
    alpha = Alpha(kk);
    fnotch = zeros(size(Lkm));
    for k= 1:length(Lkm)
        Fiber.L = Lkm(k)*1e3;
        Himdd = -abs(Fiber.Himdd(f, lambda, alpha, 'large signal')).^2;
        Himdd(f < 0) = [];
        [pks, locs] = findpeaks(Himdd);
        if isempty(locs)
            fnotch(k) = Inf;
        else
            fnotch(k) = fpos(min(locs));
        end
    end
    plot(Lkm, fnotch/1e9, 'LineWidth', 4, 'DisplayName', ['\alpha = ' num2str(alpha)])
end
plot(Lkm([1 end]), 28/2*[1 1], ':k', 'LineWidth', 4, 'DisplayName', 'Rs/2')
xlabel('Fiber length (km)', 'FontSize', 12)
ylabel('Frequency of first notch (GHz)', 'FontSize', 18)
legend('-DynamicLegend')
set(gca, 'FontSize', 18) 
% title('Location of first notch in fiber frequency response')
axis([Lkm([1 end]) 0 35])

%% Impulse response
Fiber = fiber((20/17)*1e3);
figure, hold on, box on
t = linspace(-0.1e-9, 0.1e-9, 1e3);
plot(t*1e9, real(Fiber.hdisp(t, lambda)), 'LineWidth', 4)
plot(t*1e9, -imag(Fiber.hdisp(t, lambda)), 'LineWidth', 4)
xlabel('Time (ns)', 'FontSize', 18)
ylabel('Impulse response')
legend('Filter 1', 'Filter 2')
set(gca, 'ytick', []);
set(gca, 'FontSize', 18)
grid on
