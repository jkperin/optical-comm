% Computes the coefficients of the L-point Wiener filter
% W = inv(R)*P
% The system is assumed to be of the form:
% theta(k) = theta(k-Delta) + (theta(k)-theta(k-Delta)) + n'(k)
% where phase theta evolves as a Wiener process with
% E[(theta(k)-theta(k-l))^2] = sigmapsq*|l|
% n'(k) is additive white Gaussian noise with variance sigmandsq
% The filter optimizes the estimate of theta(L-Delta) over a block of
% estimates theta(1) to theta(L)
function Wopt = wiener_fir(L, Delta, sigmandsq, sigmapsq)

% Build autocorrelation matrix of phase noise
KP = zeros(L,L);

if (Delta==L)
    r = [fliplr([1:Delta])];
else
    r = [fliplr([1:Delta]),0,[1:(L-Delta)-1]];
end;

for row=1:Delta
    for col=1:Delta
        KP(row,col) = min(r(row),r(col));
    end;
end;

for row=Delta+1:L
    for col=Delta+1:L
        KP(row,col) = min(r(row),r(col));
    end;
end;

KP = KP*sigmapsq;

% Compute autocorrelation matrix of noise
Kgamma = sigmandsq*eye(L,L);
K = KP + Kgamma;
invK = inv(K);

% Compute optimal filter coefficients
Wopt = (invK*ones(L,1))/sum(sum(invK));