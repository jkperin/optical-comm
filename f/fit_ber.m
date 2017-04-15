function [PsensdBm, f] = fit_ber(PrxdBm, BER, BERtarget)
%% Calculate receiver sensitivity i.e., receiver power PrxdBm so that BER == BERtarget
PrxdBm(BER == 0) = [];
BER(BER == 0) = [];

f = fit(PrxdBm(:), log10(BER(:)), 'linearinterp');
[PsensdBm, ~, exitflag] = fzero(@(x) f(x) - log10(BERtarget), 0);

if exitflag ~= 1
    warning('rx_sensitivity: receiver sensitivity calculation exited with exitflag = %d', exitflag);
end