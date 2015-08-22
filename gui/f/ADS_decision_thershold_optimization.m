function [Pthresh, Plevels, OptDecisionInstant] = ADS_decision_thershold_optimization(x, Pthresh, Mct, verbose)

%% Average power and extinction ratio of ADS output
n = 200*Mct+1:length(x);
nn = n(1:min(2^14, length(n)));

x = x(n);

varPl = zeros(Mct, 4);
for k = 1:Mct
    Psampled = x(k:Mct:end);

%     figure
%     hist(Psampled, 50)
%     title(sprintf('k = %d', k))
    
    Ps1 = Psampled(Psampled < Pthresh(1));
    Ps2 = Psampled(Psampled > Pthresh(1) & Psampled < Pthresh(2));
    Ps3 = Psampled(Psampled > Pthresh(2) & Psampled < Pthresh(3));
    Ps4 = Psampled(Psampled > Pthresh(3));
    
    % Must choose desicion instant that minimizes intraclass variance
    varPl(k, :) = [var(Ps1), var(Ps2), var(Ps3), var(Ps4)];
end

% [minVar, IX] = min(varPl, [], 1)
% 
% figure
% plot(varPl, '-o')
% axis([1 Mct 0 0.9])
% legend('Level 1', 'Level 2', 'Level 3', 'Level 4')

%% Chosen optimal decision instant
[~, OptDecisionInstant] = min(mean(varPl, 2));
Psampled = x(OptDecisionInstant:Mct:end);

Ps1 = Psampled(Psampled < Pthresh(1));
Ps2 = Psampled(Psampled > Pthresh(1) & Psampled < Pthresh(2));
Ps3 = Psampled(Psampled > Pthresh(2) & Psampled < Pthresh(3));
Ps4 = Psampled(Psampled > Pthresh(3));

Plevels = [mean(Ps1), mean(Ps2), mean(Ps3), mean(Ps4)];

Pthresh(1) = (Plevels(2) + Plevels(1))/2;
Pthresh(2) = (Plevels(3) + Plevels(2))/2;
Pthresh(3) = (Plevels(4) + Plevels(3))/2;


if verbose
    eyediagram(circshift(x(nn), [0, -OptDecisionInstant+1]), Mct, Mct);
    hold on
    plot((OptDecisionInstant-1)*[1 1], 15*[-1 1], 'k', 'LineWidth', 2)
    plot(floor(Mct/2)*[-1 1], Pthresh*[1 1], 'k', 'LineWidth', 2);
    axis([-8 8 0 15])
end
