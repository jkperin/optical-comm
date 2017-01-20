%% Process data saved by QPSK_BER_qsub.m
clear, clc

addpath ../
addpath ../f/
addpath ../DSP/
addpath ../../f/
addpath ../../apd/
addpath ../../soa/

folder = 'DQPSK_Analog/';

Rb = 2*112e9;
Rs = Rb/(4);
BERtarget = 1.8e-4;
Modulator = 'SiPhotonics';
ModBW = 30;
linewidth = 200;
lamb = [1250 1380];

R = 1;
q = 1.60217662e-19;

Fiber = fiber();

SNRdB2PrxdBm = @(SNRdB) 10*log10(10^(SNRdB/10)*2*q*Rs/(R*1e-3));
SNRdBref = SNRreq(BERtarget, 4, 'QAM');
PrefdBm = -35.0798; % SNRdB2PrxdBm(SNRdBref);

LineStyle = {'-', '--'};
Marker = {'o', 's', 'v'};
Color = {[51, 105, 232]/255, [153,153,155]/255, [255,127,0]/255};
Lspan = 0:0.5:10;
for l = 1:length(lamb)
    for k = 1:length(Lspan)
        filename = [folder sprintf('DQPSK_Analog_BER_L=%skm_lamb=%snm_ModBW=%sGHz_linewidth=%skHz.mat',...
        num2str(Lspan(k)), num2str(lamb(l)), num2str(ModBW), num2str(linewidth))];  
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
            [PrxdBm(l, k), ~, exitflag] = fzero(@(x) f(x) - log10(BERtarget), -28);
            figure(2), clf, hold on, box on
            hline = plot(S.Tx.PlaunchdBm, BERcount, '-o');
            plot(S.Tx.PlaunchdBm, f(S.Tx.PlaunchdBm), '-', 'Color', get(hline, 'Color'));
            axis([S.Tx.PlaunchdBm([1 end]) -8 0])
            if exitflag ~= 1
                disp('Interpolation failed')
                exitflag
                PrxdBm(l, k) = NaN;
            end
            drawnow   
            1;
        catch e
            filename
            warning(e.message)
            PrxdBm(l, k) = NaN;
        end
    end
end

D = [D(1, end:-1:1) D(2, :)];
PpendB = [PrxdBm(1, end:-1:1) PrxdBm(2, :)] -PrefdBm;

figure(1), hold on, box on
plot(D*1e6, PpendB, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'w');
xlabel('Dispersion ps/nm')
ylabel('Power penalty (dB)')
axis([-60 60 0 10])