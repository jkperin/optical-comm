%% Process data saved by QPSK_BER_qsub.m
clear, clc, close all

addpath ../f/
addpath ../DSP/
addpath ../../f/
addpath ../../apd/
addpath ../../soa/

% QPSK_BER_L=0km_SiPhotonics_BW=30GHz_DPLL_0taps_nu=200MHz_fOff=0GHz_ros=125_ENOB=Inf

BERtarget = 1.8e-4;
ros = 1.25;
nu = 200;
ENOB = [4 Inf];
fOff = 8;
CPR = {'DPLL', 'feedforward', 'feedforward'};
CPRtaps = [0, 5, 15];
Modulator = 'SiPhotonics';
ModBW = 30;
linewidth = 200;

LineStyle = {'-', '--'};
Marker = {'o', 's', 'v'};
Lspan = 0:10;
figure(1), hold on, box on
count = 1;
for enobi = 1:length(ENOB)
    for cpri = 1:length(CPR)
        Prec = zeros(size(Lspan));
        for k = 1:length(Lspan)
            filename = sprintf('QPSK_BER_L=%dkm_%s_BW=%dGHz_%s_%dtaps_nu=%dMHz_fOff=%dGHz_ros=%d_ENOB=%d',...
                Lspan(k), Modulator, ModBW, CPR{cpri}, CPRtaps(cpri), linewidth, fOff, round(100*ros), ENOB(enobi));
            
            try 
                S = load(filename);

                lber = log10(S.BER);
                valid = ~(isinf(lber) | isnan(lber));
                try
                    Prec(k) = interp1(lber(valid), S.Tx.PlaunchdBm(valid), log10(BERtarget));
                catch e
                    warning(e.message)
                    Prec(k) = NaN;
                    lber(valid)
                end 
            catch e
                warning(e.message)
                Prec(k) = NaN;
                filename
            end
            
            figure(2), box on, hold on
            plot(S.Tx.PlaunchdBm, lber, '-o')
        end
        
        figure(1)
        hline(count) = plot(Lspan, Prec, 'LineStyle', LineStyle{enobi}, 'Marker', Marker{cpri});
        count = count + 1;
        drawnow
    end
end

xlabel('Fiber length (km)')
ylabel('Receiver sensitivity (dBm)')