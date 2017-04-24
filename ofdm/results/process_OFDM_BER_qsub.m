%% Process data saved by OFDM_BER_qsub.m
clear, clc, close all

addpath ../
addpath ../../f/
addpath ../../apd/

folder = '';

BERtarget = 1.8e-4;
OFDMtype = {'DC-OFDM'};
Amplified = {0};
ModBWGHz = 30;
ENOB = [5];

LineStyle = {'-', '--'};
Marker = {'o', 's', 'v'};
Color = {[51, 105, 232]/255, [153,153,155]/255, [255,127,0]/255};
Lkm = 0:1:15;

Fiber = fiber();
D = zeros(1, length(Lkm));
PrxdBm = zeros(length(OFDMtype), length(Amplified), length(Lkm));
for n = 1:length(OFDMtype)
    for a = 1:length(Amplified)
        for k = 1:length(Lkm)
            filename = [folder sprintf('%s_BER_Amplified=%d_L=%dkm_ModBW=%dGHz_ENOB=%d.mat',...
                OFDMtype{n}, Amplified{a}, Lkm(k), ModBWGHz, ENOB(n))];  
            try 
                S = load(filename, '-mat');
                D(k) = Fiber.D(S.Tx.Laser.wavelength)*S.Fiber.L/1e3;

                BERcount = 0;
                BERtheory = 0;
                BERgauss = 0;
                counter = 0;
                for i = 1:S.sim.Realizations
                    BERcount = BERcount + S.BER(i).count;
                    BERtheory = BERtheory + S.BER(i).theory;
                    BERgauss = BERgauss + S.BER(i).gauss;
                    counter = counter + 1;
                end
                BERcount = BERcount/counter;
                BERtheory = BERtheory/counter;
                BERgauss = BERgauss/counter;

                BERcountlog = log10(BERcount);
                
                
                idx = find(BERcountlog <= -2 & BERcountlog >= -5.5);
                [PrxdBm(n, a, k), f] = fit_ber(S.Tx.PtxdBm(idx), BERcount(idx), BERtarget);
                
                % Plot
                figure(2), clf, hold on, box on
                hline = plot(S.Tx.PtxdBm, BERcountlog, '-o');
                plot(S.Tx.PtxdBm, f(S.Tx.PtxdBm), '-', 'Color', get(hline, 'Color'));
                plot(S.Tx.PtxdBm, log10(BERgauss), '--', 'Color', get(hline, 'Color'));
                plot(S.Tx.PtxdBm, log10(BERtheory), ':k');
                axis([S.Tx.PtxdBm([1 end]) -8 0])
                title(sprintf('L = %.1f km, D = %.2f', S.Fiber.L/1e3, D(k)*1e6))
                drawnow   
            catch e
                filename
                warning(e.message)
                PrxdBm(n, a, k) = NaN;
            end
        end
    end
    figure(1), hold on, box on
    plot(D*1e6, squeeze(PrxdBm(n, a, :)), '-or', 'LineWidth', 2,...
        'MarkerFaceColor', 'w', 'DisplayName', sprintf('%s, %d', OFDMtype{n}, Amplified{a}));
end

xlabel('Dispersion (ps/nm)')
ylabel('Receiver sensitivity (dBm)')

