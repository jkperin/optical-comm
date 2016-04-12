%% Process data saved by QPSK_BER_qsub.m
clear, clc, close all

addpath ../f/
addpath ../DSP/
addpath ../../f/
addpath ../../apd/
addpath ../../soa/

% QPSK_BER_L=0km_SiPhotonics_BW=30GHz_nu=200MHz_fOff=0GHz_ros=125_ENOB=Inf

BERtarget = 1.8e-4;
ros = 1.25;
nu = 200;
ENOB = [4 Inf];
fOff = 0;
Modulator = 'SiPhotonics';
ModBW = 40;
linewidth = 200;

LineStyle = {'-', '--'};
Marker = {'o', 's', 'v'};
Color = {[51, 105, 232]/255, [255,127,0]/255};
Lspan = 0:10;
figure(1), hold on, box on
count = 1;
for enobi = 1:length(ENOB)
    Prec = zeros(size(Lspan));
    for k = 1:length(Lspan)
        filename = sprintf('DQPSK_BER_L=%dkm_%s_BW=%dGHz_nu=%dkHz_fOff=%dGHz_ros=%d_ENOB=%d',...
            Lspan(k), Modulator, ModBW, linewidth, fOff, round(100*ros), ENOB(enobi));

        try 
            S = load(filename);

            lber = log10(S.BER.count);
            valid = ~(isinf(lber) | isnan(lber)) & lber > -4.5 & lber < -3;
            try
                Prec(k) = spline(lber(valid), S.Tx.PlaunchdBm(valid), log10(BERtarget));
            catch e
                warning(e.message)
                Prec(k) = NaN;
                lber(valid)
            end 

            lber_theory = log10(S.BER.theory);
            Prec_theory(k) = spline(lber_theory, S.Tx.PlaunchdBm, log10(BERtarget));

            figure(2), box on, hold on
            h = plot(S.Tx.PlaunchdBm, lber, '-o');
            plot(Prec(k), log10(BERtarget), '*', 'Color', get(h, 'Color'), 'MarkerSize', 6)
            plot(S.Tx.PlaunchdBm, lber_theory, 'k')
            axis([S.Tx.PlaunchdBm([1 end]) -8 0])
        catch e
            warning(e.message)
            Prec(k) = NaN;
            Prec_theory(k) = NaN;
            filename
        end
    end
    Prec_dqpsk_count{enobi} = Prec;
    Prec_dqpsk_theory = Prec_theory;
    
    figure(1)
    hline(count) = plot(Lspan, Prec, 'Color', Color{1}, 'LineStyle', LineStyle{enobi}, 'LineWidth', 2);
    plot(Lspan, Prec_theory, 'k', 'LineWidth', 2)
    count = count + 1;
    drawnow
end

xlabel('Fiber length (km)')
ylabel('Receiver sensitivity (dBm)')
% axis([0 10 -33 -27])