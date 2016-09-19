function [lambchs, Disp] = optimal_wdm_channels(Dlamb, Nch, verbose)
%% Find set of Nch WDM channels with miminum sum of dispersion
% Inputs:
% - Dlamb: wavelength spacing in m
% - Nch: number of WDM channels
% Outputs:
% - lambchs: wavelengths of the WDM channels
% - Disp: corresponding dispersion at each wavelength in s/m^2

addpath ../../f % class fiber path

Fiber = fiber();

chs = 0:Nch-1;

[x0,~,exitflag] = fminbnd(@(x0) max(abs(Fiber.D(x0*1e-9 + Dlamb*chs))), 1200, 1500);

if exitflag ~= 1
    warning(sprintf('optimal_wdm_channels/WDM channel optimization ended with exit flag %d\n', exitflag));
end

lambchs = x0*1e-9 + Dlamb*chs;
Disp = Fiber.D(lambchs);

% Plot
if exist('verbose', 'var') && verbose
    lamb = linspace(lambchs(1), lambchs(end));
    figure, hold on, box on
    plot(lamb*1e9, 1e6*Fiber.D(lamb), 'k', 'LineWidth', 2)
    stem(lambchs*1e9, 1e6*Disp, 'k', 'LineWidth', 2, 'MarkerFaceColor', 'w')
    xlabel('Wavelength (nm)', 'FontSize', 12)
    ylabel('Dispersion (ps/(nm\cdot km))', 'FontSize', 12)
    set(gca, 'FontSize', 12)
end