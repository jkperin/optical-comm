%% Number of taps required by equalization
clear, clc, close all

addpath ../../f/
addpath ../f/

[lambchs20, Disp20] = optimal_wdm_channels(20e-9, 8);
[~, idx] = max(abs(Disp20));
wavelengts(1) = lambchs20(idx);

[lambchs45, Disp45] = optimal_wdm_channels(4.5e-9, 8);
[~, idx] = max(abs(Disp45));
wavelengts(2) = lambchs45(idx);

wavelengts(3) = 1550e-9;

ros = 5/4;
Rs = 56e9;

Lspan = 0:0.05:30;

Fiber = fiber();
Fiber.PMD = true;

for l = 1:length(Lspan)
    Fiber.L = Lspan(l)*1e3;
    
    for k = 1:length(wavelengts)
        [Ncd(k, l), Npmd(k, l)] = Fiber.Ntaps(Rs, ros, wavelengts(k));
    end
end
    
figure, hold on, box on
plot(Lspan, ceil(Ncd(1, :)), 'LineWidth', 2)
plot(Lspan, ceil(Ncd(2, :)), 'LineWidth', 2)
plot(Lspan, ceil(Ncd(3, :)), 'LineWidth', 2)
xlabel('Fiber length (km)', 'FontSize', 12)
ylabel('Required number of taps', 'FontSize', 12)
leg = legend('\Delta\lambda = 20 nm: \lambda = 1390 nm', '\Delta\lambda = 4.5 nm: \lambda = 1328 nm', '\lambda = 1550 nm')
leg.FontSize = 12;
set(gca, 'FontSize', 12)
grid on