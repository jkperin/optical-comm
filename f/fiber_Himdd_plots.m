%% Fiber IM-DD frequency response
clear, clc, close all

f = linspace(0, 112, 1e3)*1e9;
lambda = 1380e-9;

Fiber = fiber(10e3);
%% Large signal, fixed length, variable chirp
figure, hold on, box on
plot(f/1e9, abs(Fiber.Himdd(f, lambda, 0, 'large signal')).^2, 'LineWidth', 2)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -1, 'large signal')).^2, 'LineWidth', 2)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -2, 'large signal')).^2, 'LineWidth', 2)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -3, 'large signal')).^2, 'LineWidth', 2)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -4, 'large signal')).^2, 'LineWidth', 2)
xlabel('Frequency (GHz)', 'FontSize', 12)
ylabel('|H_{IM-DD}(f)|^2', 'FontSize', 12)
leg = legend('\alpha = 0', '-1', '-2', '-3', '-4');
set(leg, 'FontSize', 12);
set(gca, 'FontSize', 12)
% title(sprintf('Large-signal IM-DD frequency response for L = %.2f km', Fiber.L/1e3))
% saveas(gca, sprintf('Hfiber_large_signal_L=%dkm.png', Fiber.L/1e3))

%% Small signal, fixed length, variable chirp
figure, hold on, box on
plot(f/1e9, abs(Fiber.Himdd(f, lambda, 0, 'small signal')).^2, 'LineWidth', 2)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -1, 'small signal')).^2, 'LineWidth', 2)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -2, 'small signal')).^2, 'LineWidth', 2)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -3, 'small signal')).^2, 'LineWidth', 2)
plot(f/1e9, abs(Fiber.Himdd(f, lambda, -4, 'small signal')).^2, 'LineWidth', 2)
xlabel('Frequency (GHz)', 'FontSize', 12)
ylabel('|H_{IM-DD}(f)|^2', 'FontSize', 12)
leg = legend('\alpha = 0', '-1', '-2', '-3', '-4');
set(leg, 'FontSize', 12);
set(gca, 'FontSize', 12)
% title(sprintf('Small-signal IM-DD frequency response for L = %.2f km', Fiber.L/1e3))
% saveas(gca, sprintf('Hfiber_small_signal_L=%dkm.png', Fiber.L/1e3))

%% Large signal, variable length, fixed chirp
figure, hold on, box on
alpha = 0;
Fiber.L = 1e3;
plot(f/1e9, abs(Fiber.Himdd(f, lambda, alpha, 'large signal')).^2, 'LineWidth', 2, 'DisplayName', sprintf('%.2f ps/nm', 1e6*abs(Fiber.D(lambda))*Fiber.L/1e3))
Fiber.L = 5e3;
plot(f/1e9, abs(Fiber.Himdd(f, lambda, alpha, 'large signal')).^2, 'LineWidth', 2, 'DisplayName', sprintf('%.2f ps/nm', 1e6*abs(Fiber.D(lambda))*Fiber.L/1e3))
Fiber.L = 10e3;
plot(f/1e9, abs(Fiber.Himdd(f, lambda, alpha, 'large signal')).^2, 'LineWidth', 2, 'DisplayName', sprintf('%.2f ps/nm', 1e6*abs(Fiber.D(lambda))*Fiber.L/1e3))
Fiber.L = 15e3;
plot(f/1e9, abs(Fiber.Himdd(f, lambda, alpha, 'large signal')).^2, 'LineWidth', 2, 'DisplayName', sprintf('%.2f ps/nm', 1e6*abs(Fiber.D(lambda))*Fiber.L/1e3))
Fiber.L = 20e3;
plot(f/1e9, abs(Fiber.Himdd(f, lambda, alpha, 'large signal')).^2, 'LineWidth', 2, 'DisplayName', sprintf('%.2f ps/nm', 1e6*abs(Fiber.D(lambda))*Fiber.L/1e3))
xlabel('Frequency (GHz)', 'FontSize', 12)
ylabel('|H_{IM-DD}(f)|^2', 'FontSize', 12)
legend('-DynamicLegend')
set(leg, 'FontSize', 12);
set(gca, 'FontSize', 12)
% title(sprintf('Large-signal IM-DD frequency response for alpha = %.1f km', alpha))
% saveas(gca, sprintf('Hfiber_large_signal_alpha=%dkm.png', alpha))

%% First notch location
Lkm = 0:50;
Alpha = 0:-1:-3;
fpos = f(f >= 0);
figure, hold on, box on
D = zeros(size(Lkm));
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
        D(k) = 1e6*Fiber.D(lambda)*Fiber.L/1e3;
    end
    leg(k) = plot(abs(D), fnotch/1e9, 'LineWidth', 2, 'DisplayName', ['\alpha = ' num2str(alpha)]);
end
plot(abs(D([1 end])), 112/2*[1 1], ':k', 'DisplayName', '56 GHz','LineWidth', 2)
plot(abs(D([1 end])), 56/2*[1 1], ':k', 'DisplayName', '28 GHz', 'LineWidth', 2)
xlabel('Dispersion (ps/nm)', 'FontSize', 12)
ylabel('First notch frequency (GHz)', 'FontSize', 12)
legend('-DynamicLegend')
set(gca, 'FontSize', 12) 
% title('Location of first notch in fiber frequency response')
% axis([0 180 20 80])
m = matlab2tikz(gca);
m.write('first-notch-freq.tex')

%% Impulse response
Fiber = fiber((20/17)*1e3);
figure, hold on, box on
t = linspace(-0.1e-9, 0.1e-9, 1e3);
plot(t*1e9, real(Fiber.hdisp(t, lambda)), 'LineWidth', 2)
plot(t*1e9, -imag(Fiber.hdisp(t, lambda)), 'LineWidth', 2)
xlabel('Time (ns)', 'FontSize', 12)
ylabel('Impulse response')
legend('Filter I', 'Filter Q')
set(gca, 'ytick', []);
set(gca, 'FontSize', 12)
a = axis;
axis([-0.02 0.02 a(3:4)]);
grid on