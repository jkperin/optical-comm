function [xaligned, lag] = align_signals(x, xref)
%% Delay/Advance signal x in order to maximizer cross-correlation between two signals
[c, lags] = xcorr(xref, x, 20);

[~, ind] = max(abs(c));

lag = lags(ind)
xaligned = circshift(x, [0 lag]);
1;

if abs(lag) > 3
    figure, plot(lags, abs(c))
    1;
end