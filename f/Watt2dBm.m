function PdBm = Watt2dBm(P)
% Converts power from Watt to dBm
PdBm = 10*log10(P/1e-3);