function Ysq = Quant(Y,Quantization)

    Y1i = real(Y(1,:));
    Y1q = imag(Y(1,:));
    Y2i = real(Y(2,:));
    Y2q = imag(Y(2,:));

    % quantize
    Bits = Quantization.Bits;
    Range = Quantization.Range;
    Y1iq = quantz(Y1i,Range,Bits);   % quantized samples of inphase in x pol.
    Y1qq = quantz(Y1q,Range,Bits);   % quantized samples of quad in x pol.
    Y2iq = quantz(Y2i,Range,Bits);   % quantized samples of inphase in y pol.
    Y2qq = quantz(Y2q,Range,Bits);   % quantized samples of quad in y pol.

    Ysq = [Y1iq + 1i*Y1qq; Y2iq + 1i*Y2qq];
end
function quantout = quantz(quantin, quantrange, quantbits)
    quantlvls = 2^quantbits;
    deltaquant = quantrange/(quantlvls-1);
    quantout = max(-quantlvls/2*deltaquant,min((quantlvls/2-1)*deltaquant,deltaquant*round(quantin/deltaquant)));
end