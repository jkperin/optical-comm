function SNRdB = SNRreq(BER, M, format)
%% SNR required to achieve BER using uncoded M-{qam, dpsk, etc} modulation in AWGN channel
% Note that this function assumes the relation: SNR = log2(M)Eb/N0, which 
% is valid for complex signals since it assumes that the  noise power is
% 2N0 i.e., two components whose two-sided noise PSD is N0/2. For a real
% signal, the relation should be SNR = 2log2(M)Eb/N0.

[SNRdB, ~, exitflag] = fzero(@(X) log10(berawgn(X - 10*log10(log2(M)), lower(format), M)) - log10(BER), 10);

if exitflag ~= 1
    warning('SNRreq: could not find SNR to achieve target BER')
    SNRdB = NaN;
end