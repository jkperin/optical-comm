function [Pthresh, Plevels, OptDecisionInstant] = ADS_decision_thershold_optimization(x, mpam, Mct, verbose)

%% Average power and extinction ratio of ADS output
Pavg = 8; % average power in mW
rexdB = -5; % extinction ratio in dB

n = 100*Mct+1:length(x)-100*Mct;
nn = n(1:min(2^14, length(n)));

x = x(n);

[~, Pthresh] = mpam.adjust_levels(Pavg, rexdB);

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
    eyediagram(x(nn), Mct, Mct);
    hold on
    plot((OptDecisionInstant-1)*[1 1], 15*[-1 1], 'k', 'LineWidth', 2)
    plot(ceil([-Mct Mct]/2), Pthresh*[1 1], 'k', 'LineWidth', 2);
end