%% Process data saved by QPSK_BER_qsub.m
clear, clc

addpath ../
addpath ../f/
addpath ../DSP/
addpath ../../f/
addpath ../../apd/
addpath ../../soa/

% QPSK_Analog_BER_L=0.5km_lamb=1380nm_ModBW=30GHz_OPLL-costas_Npol=1_linewidth=200kHz_delay=250ps

Rb = 2*112e9;
Rs = Rb/(4);
BERtarget = 1.8e-4;
CPR = {'costas', 'logic'};
delay = 250;
Modulator = 'SiPhotonics';
ModBW = 30;
linewidth = 200;
lamb = [1250 1380];

R = 1;
q = 1.60217662e-19;

SNRdB2PrxdBm = @(SNRdB) 10*log10(10^(SNRdB/10)*2*q*Rs/(R*1e-3));
SNRdBref = SNRreq(BERtarget, 4, 'QAM');
PrefdBm = -35.0798; % SNRdB2PrxdBm(SNRdBref);

LineStyle = {'-', '--'};
Marker = {'o', 's', 'v'};
Color = {[51, 105, 232]/255, [153,153,155]/255, [255,127,0]/255};
Lspan = 0:0.5:10;

Fiber = fiber();
PrxdBm = zeros(2, 1, length(lamb), length(Lspan));
D = zeros(length(lamb), length(Lspan));
for n = 1:2
    for Ncpr = 1
        for l = 1:length(lamb)
            for k = 1:length(Lspan)
                filename = sprintf('QPSK_Analog_BER_L=%skm_lamb=%snm_ModBW=%sGHz_OPLL-%s_Npol=%s_linewidth=%skHz_delay=%sps.mat',...
                num2str(Lspan(k)), num2str(lamb(l)), num2str(ModBW), CPR{n}, num2str(Ncpr), num2str(linewidth),...
                num2str(delay));  
                try 
                    S = load(filename, '-mat');
                    D(l, k) = Fiber.D(S.Tx.Laser.wavelength)*S.Fiber.L/1e3;

                    BERcount = 0;
                    BERtheory = 0;
                    counter = 0;
                    for i = 1:S.sim.Realizations
                        BERcount = BERcount + S.BER(i).count;
                        BERtheory = BERtheory + S.BER(i).theory;
                        counter = counter + 1;
                    end
                    BERcount = log10(BERcount/counter);
                    BERtheory = log10(BERtheory/counter);

                    idx = find(BERcount <= -2 & BERcount >= -5.5);
                    f = fit(S.Tx.PlaunchdBm(idx).', BERcount(idx).', 'linearinterp');
                    [PrxdBm(n, Ncpr, l, k), ~, exitflag] = fzero(@(x) f(x) - log10(BERtarget), -28);
                    figure(2), clf, hold on, box on
                    hline = plot(S.Tx.PlaunchdBm, BERcount, '-o');
                    plot(S.Tx.PlaunchdBm, f(S.Tx.PlaunchdBm), '-', 'Color', get(hline, 'Color'));
                    axis([S.Tx.PlaunchdBm([1 end]) -8 0])
                    if exitflag ~= 1
                        disp('Interpolation failed')
                        exitflag
                        PrxdBm(n, Ncpr, l, k) = NaN;
                    end
                    drawnow   
                    1;
                catch e
                    filename
                    warning(e.message)
                    PrxdBm(n, Ncpr, l, k) = NaN;
                end
            end
        end
        Dcomb = [D(1, end:-1:1) D(2, :)];
        PpendBcomb(n, Ncpr, :) = [squeeze(PrxdBm(n, Ncpr, 1, end:-1:1)).' squeeze(PrxdBm(n, Ncpr, 2, :)).'] -PrefdBm;
        figure(1), hold on, box on
        plot(Dcomb*1e6, squeeze(PpendBcomb(n, Ncpr, :)), '-or', 'LineWidth', 2, 'MarkerFaceColor', 'w');
    end
end

xlabel('Dispersion ps/nm')
ylabel('Receiver sensitivity (dBm)')
axis([-60 60 0 10])

