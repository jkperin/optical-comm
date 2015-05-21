function quantout = quantz(quantin, quantrange, numquantbits);
numquantlvls = 2^numquantbits;
deltaquant = quantrange/(numquantlvls-1);
quantout = max(-numquantlvls/2*deltaquant,min((numquantlvls/2-1)*deltaquant,deltaquant*round(quantin/deltaquant)));
