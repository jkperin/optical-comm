%% Process data saved by QPSK_BER_qsub.m
clear, clc, close all

addpath ../f/
addpath ../DSP/
addpath ../../f/
addpath ../../apd/
addpath ../../soa/

BERtarget = 1.8e-4;
Modulator = 'SiPhotonics';
ModBW = 40;
linewidth = 200;
ideal= [1 0 0 0 0];
Delay = [0 0 40 60 80];

LineStyle = {'-', '--'};
Marker = {'o', 's', 'v'};
Color = {[51, 105, 232]/255, [228,26,28]/255, [0,153,37]/255, [255,127,0]/255,...
            [152,78,163]/255, [166,86,40]/255, [153,153,155]/255, [77,175,74]/255, [0,0,0]};  
Lspan = 0:10;
figure(1), hold on, box on
count = 1;
for i = 1:length(ideal)
    Prec = zeros(size(Lspan));
    for k = 1:length(Lspan)
         filename = sprintf('QPSK_Analog_BER_L=%dkm_%s_BW=%dGHz_linewidth=%dkHz_ideal=%d_delay=%dsamples',...
            Lspan(k), Modulator, ModBW, linewidth, ideal(i), Delay(i));

        try 
            S = load(filename);

            lber = log10(S.BER.count);
            valid = ~(isinf(lber) | isnan(lber)) & lber > -4.5;
            try
                [~, idx] = unique(lber(valid));
                f = fit(lber(valid(idx)).', S.Tx.PlaunchdBm(valid(idx)).', 'linearinterp');
                Prec(k) = f(log10(BERtarget));
            catch e
                warning(e.message)
                Prec(k) = NaN;
                lber(valid)
            end 

            lber_theory = log10(S.BER.theory(valid));
            Prec_theory(k) = interp1(lber_theory, S.Tx.PlaunchdBm(valid), log10(BERtarget));

            figure(2), box on, hold on
            h = plot(S.Tx.PlaunchdBm, lber, '-o');
            plot(Prec(k), log10(BERtarget), '*', 'Color', get(h, 'Color'), 'MarkerSize', 6)
            plot(S.Tx.PlaunchdBm(valid), lber_theory, 'k')
            axis([S.Tx.PlaunchdBm([1 end]) -8 0])
        catch e
            warning(e.message)
            Prec(k) = NaN;
            Prec_theory(k) = NaN;
            filename
        end
    end
    Prec_qpsk_theory =  Prec_theory;
    Prec_qpsk_count{i} = Prec;

    figure(1), hold on
    hline(count) = plot(Lspan, Prec, 'Color', Color{i}, 'LineStyle', LineStyle{1}, 'Marker', Marker{1}, 'LineWidth', 2);
    count = count + 1;
    drawnow
end

hline(end+1) = plot(Lspan, Prec_theory, 'k', 'LineWidth', 2)
xlabel('Fiber length (km)')
ylabel('Receiver sensitivity (dBm)')
legend(hline, 'Ideal components', 'Loop delay = 30ps', 'Loop delay = 90ps', 'AWGN limit') 
% axis([0 10 -33 -27])
% legend(hline, 'DPLL', 'FIR Feedforward w/ 15 taps')