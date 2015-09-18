function [Pthresh, Plevels, OptDecisionInstant] = PAM4_timing_opt(x, Mct, verbose)
%% Calculate decision thresholds, levels, and optimal decision instant for 4-PAM signal x
% Given x the levels are determined by calculating the histogram and
% finding the 4 highest peaks in the histogram, which corresponds to the 4
% PAM levels. The decision thresholds are calculated as the average of the
% PAM levels, because equally spaced levels are assumed.
% The optimal decision instant is calculated by selecting the instant point
% that minimizes the variance of each level i.e., the error probability
M = 4;

%% Find levels
[c, bins] = hist(x, 100);
bint = linspace(bins(1), bins(end), 1e3);
cc = spline(bins, c, bint);

[pks, locs] = findpeaks(cc);

[~, ix] = sort(pks, 'descend');

Plevels = bint(locs(ix(1:M)));

Plevels = sort(Plevels);

if verbose
    figure, hold on
    bar(bins, c)
    plot(bint, cc)
    plot(bint(locs), cc(locs), '*g')
    plot(bint(locs(ix(1:M))), cc(locs(ix(1:M))), 'sr', 'MarkerSize', 6)
end

%% Calculate Decision thresholds
Pthresh(1) = mean(Plevels([1 2]));
Pthresh(2) = mean(Plevels([2 3]));
Pthresh(3) = mean(Plevels([3 4]));
    
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

%% Choose optimal decision instant
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
