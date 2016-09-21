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
CPRtype = {'Costas', 'Logic'};
linewidth = 200;
ideal= [1 1 1];
Delay = [0 250 500];

LineStyle = {'-', '--'};
Marker = {'o', 'o', 'v'};
Color = {[51, 105, 232]/255, [228,26,28]/255, [0,153,37]/255, [255,127,0]/255,...
            [152,78,163]/255, [166,86,40]/255, [153,153,155]/255, [77,175,74]/255, [0,0,0]};  
Lspan = 0:10;
figure(1), hold on, box on
count = 1;
for ind1 = 1:length(CPRtype)
    for i = 1:length(ideal)
        Prec = zeros(size(Lspan));
        for k = 1:length(Lspan)
             filename = sprintf('analog/Analog_QPSK_BER_%s_%s_L=%dkm_%s_BW=%dGHz_linewidth=%dkHz_ideal=%d_delay=%dps',...
                'EPLL', CPRtype{ind1}, Lspan(k), Modulator, ModBW, linewidth, ideal(i), Delay(i));

            try 
                S = load(filename);

                lber = log10(S.BER.count);
                valid = ~(isinf(lber) | isnan(lber)) & lber > -4.5;
                try
                    Pval = S.Tx.PlaunchdBm(valid).';
                    lber_valid = lber(valid);
                    [~, ia] = unique(lber_valid);
                    lber_unique = lber_valid(ia).';
                    Punique = Pval(ia);
                    f = fit(lber_unique, Punique, 'linearinterp');
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
        hline(count) = plot(Lspan, Prec, 'Color', Color{i}, 'LineStyle', LineStyle{ind1}, 'Marker', Marker{ind1}, 'LineWidth', 2, 'MarkerFaceColor', 'w');
        count = count + 1;
        drawnow
    end
end

hline(end+1) = plot(Lspan, Prec_theory, 'k', 'LineWidth', 2);
xlabel('Fiber length (km)', 'FontSize', 12)
ylabel('Receiver sensitivity (dBm)', 'FontSize', 12)
set(gca, 'FontSize', 12)

axis([0 10 -34 -24])
% axis([0 10 -33 -27])
% legend(hline, 'DPLL', 'FIR Feedforward w/ 15 taps')

Preq_QPSK_eq3taps_4ENOB = [-30.4164745074869,-30.4058273271944,-30.3330441704541,-30.3715741109382,-30.2602322104346,-30.2036180718854,-29.9936936068042,-29.9671004192639,-29.7529982889092,-29.4985634102842,-29.3006527908808];
Preq_QPSK_eq7taps_4ENOB = [-31.4125107114068,-31.3779270998451,-31.3164488824664,-31.2780039701913,-31.3520527618346,-31.2563675630953,-31.0883355328784,-31.2731071351150,-31.1648786133992,-31.0440648321642, NaN];  
Preq_QPSK_theory = [-33.0063808502409,-33.0063808502409,-33.0063808502409,-33.0063808502409,-33.0063808502409,-33.0063808502409,-33.0063808502409,-33.0063808502409,-33.0063808502409,-33.0063808502409,-33.0063808502409];

h1 = plot(Lspan, Preq_QPSK_eq3taps_4ENOB, 'Color', [153,153,155]/255, 'LineStyle', '-', 'Marker', 's', 'LineWidth', 2, 'MarkerFaceColor', 'w');
h2 = plot(Lspan, Preq_QPSK_eq7taps_4ENOB, 'Color', [153,153,155]/255, 'LineStyle', '--', 'Marker', 's', 'LineWidth', 2, 'MarkerFaceColor', 'w');
h3 = plot(Lspan, Preq_QPSK_theory, 'k', 'LineWidth', 2)
hline = [h1 h2 hline(1:3) h3];
leg = legend(hline, 'DSP w/ 3-tap equalizer', 'DSP w/ 7-tap equalizer', 'Loop delay = 0 ps', '250 ps', '500 ps', 'AWGN limit', 'Location', 'NorthWest');
set(leg, 'FontSize', 12)