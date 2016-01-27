function sigmaesq = phase_error_fir(L,Delta,sigmandsq,sigmapsq,W)

% Compute theoretical noise performance
sigmaesq = sigmandsq*sum(W.^2);

for m = 0:Delta-1
    sigmaesq = sigmaesq + sigmapsq*(sum(W(1:m+1))^2);
end;

for m = Delta:L-2
    sigmaesq = sigmaesq + sigmapsq*(sum(W(m+2:L))^2);
end;