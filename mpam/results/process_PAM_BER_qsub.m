%% Process data saved by PAM_BER_qsub.m
clear, clc, close all

addpath ../
addpath ../../f/
addpath ../../apd/

folder = '';

% PAM_BER_L=0.5km_lamb=1380nm_ModBW=30GHz_amplified=0_Ntaps=9_ENOB=5_ros=1.25

BERtarget = 1.8e-4;
Amplified = {0, 1};
ModBWGHz = 30;
ENOB = 5;
Ntaps = 9;
lamb = 1380;
ros = 1.25;

LineStyle = {'-', '--'};
Marker = {'o', 's', 'v'};
Color = {[51, 105, 232]/255, [153,153,155]/255, [255,127,0]/255};
Lkm = 0:0.5:10;

Fiber = fiber();
D = zeros(1, length(Lkm));
PrxdBm = zeros(length(Amplified), length(Lkm));
for a = 1:length(Amplified)
    for k = 1:length(Lkm)
        filename = [folder sprintf('PAM_BER_L=%skm_lamb=%dnm_ModBW=%dGHz_amplified=%d_Ntaps=%d_ENOB=%d_ros=%s.mat',...
            num2str(Lkm(k)), lamb, ModBWGHz, Amplified{a}, Ntaps, ENOB, num2str(ros))];  
        try 
            S = load(filename, '-mat');
            D(k) = 1e6*Fiber.D(S.Tx.Laser.wavelength)*S.Fiber.L/1e3;

            % Realizations were already averaged
            BERcount = log10(S.BER.count);
            BERgauss = log10(S.BER.gauss);

            idx = find(BERcount <= -2 & BERcount >= -5.5);
            f = fit(S.Tx.PtxdBm(idx).', BERcount(idx).', 'linearinterp');
            [PrxdBm(a, k), ~, exitflag] = fzero(@(x) f(x) - log10(BERtarget), -28);
            figure(2), clf, hold on, box on
            hline = plot(S.Tx.PtxdBm, BERcount, '-o');
            plot(S.Tx.PtxdBm, f(S.Tx.PtxdBm), '-', 'Color', get(hline, 'Color'));
            plot(S.Tx.PtxdBm, BERgauss, '--', 'Color', get(hline, 'Color'));
            axis([S.Tx.PtxdBm([1 end]) -8 0])
            title(sprintf('L = %.1f km, D = %.2f', S.Fiber.L/1e3, D(k)))
            if exitflag ~= 1
                disp('Interpolation failed')
                exitflag
                PrxdBm(n, a, k) = NaN;
            end
            drawnow   
            1;
        catch e
            filename
            warning(e.message)
            PrxdBm(a, k) = NaN;
        end
    end
    figure(1), hold on, box on
    plot(D, PrxdBm(a, :), '-or', 'LineWidth', 2,...
        'MarkerFaceColor', 'w', 'DisplayName', sprintf('Amplified = %d', Amplified{a}));
    xlabel('Dispersion (ps/nm)')
    ylabel('Receiver sensitivity (dBm)')
end

