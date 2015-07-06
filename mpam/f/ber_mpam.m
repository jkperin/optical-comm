%% Calculate BER of M-PAM signal with signal-dependent or signal-independent AWGN
% Plevels = amplitude of each level
% Pthresh = decision thresholds
% var_levels = noise variance assuming each level was transmitted
function ber = ber_mpam(Plevels, Pthresh, var_levels)

M = length(Plevels);
ser = 0;
for k = 1:M
    if k == 1
        ser = ser + 1/M*qfunc((Pthresh(1) - Plevels(1))/sqrt(var_levels(1)));
    elseif k == M
        ser = ser + 1/M*qfunc((Plevels(M) - Pthresh(M-1))/sqrt(var_levels(M)));
    else
        ser = ser + 1/M*qfunc((Pthresh(k) - Plevels(k))/sqrt(var_levels(k)));
        ser = ser + 1/M*qfunc((Plevels(k) - Pthresh(k-1))/sqrt(var_levels(k)));
    end
end

ber = ser/log2(M);
