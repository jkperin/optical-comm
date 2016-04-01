function SNRdB = SNRreq(BER, M, format)
%% SNR required to achieve BER using uncoded M-{qam, dpsk, etc} modulation in AWGN channel

[SNRdB, ~, exitflag] = fzero(@(X) log10(berawgn(X - 10*log10(log2(M)), lower(format), M)) - log10(BER), 10);

if exitflag ~= 1
    error('SNRreq/could not find SNR to achieve target BER')
end