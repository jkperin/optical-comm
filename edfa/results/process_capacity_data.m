%% 
clear, clc, close all

addpath ../

edf_type = 'principles_type2';
totalLengthKm=14e3;
totalPump = 5:10;
spanLengthKm = 30:5:80;

for p = 1:length(totalPump)
    for l = 1:length(spanLengthKm)
        try
            S = load(sprintf('principles_type2/capacity_%s_totalPump=%dW_Ltotal=%d_Lspan=%dkm.mat',...
            edf_type, totalPump(p), totalLengthKm, spanLengthKm(l)));

            Lopt(p, l) = S.E.L;
            C(p, l) = S.C;
        catch e
            disp(e.message)
            Lopt(p, l) = NaN;
            C(p, l) = NaN;
        end
    end
    
    figure(1), hold on, box on
    plot(spanLengthKm, Lopt(p, :), 'DisplayName', sprintf('Total pump = %d', totalPump(p)))
    xlabel('Span length (km)')
    ylabel('Optimal EDF length (m)')
    legend('-dynamiclegend')
       
    figure(2), hold on, box on
    plot(spanLengthKm, C(p, :), 'DisplayName', sprintf('Total pump = %d', totalPump(p)))
    xlabel('Span length (km)')
    ylabel('Total spectrum efficiency (bits/s/Hz)')
    legend('-dynamiclegend')
    axis([30 80 50 140])
    
end
