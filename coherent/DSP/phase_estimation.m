function [Y, phi] = phase_estimation(X, eq, Ytrain, sim)  

M = sim.M;
Ntrain = eq.Ntrain;
mu = eq.mu;
Nsymb = length(X);

% Initialize filters coefficients
phix = zeros(1, Nsymb+1);
phiy = zeros(1, Nsymb+1);

yx_hat = zeros(1, Nsymb); % output in x pol stream
yy_hat = zeros(1, Nsymb); % output in y pol stream
for n = 1:Nsymb % n runs over symbols
    % Phase tracker
    yx_hat(n) = X(1, n)*exp(-1j*phix(n));
    yy_hat(n) = X(2, n)*exp(-1j*phiy(n));

    if n < Ntrain % Training based on training sequence
        ephix = yx_hat(n) - Ytrain(1, n);
        ephiy = yy_hat(n) - Ytrain(2, n);
    else % based on deteced symbols
        ephix = yx_hat(n) - qammod(qamdemod(yx_hat(n), M, 0, 'Gray'), M, 0, 'Gray');
        ephiy = yy_hat(n) - qammod(qamdemod(yy_hat(n), M, 0, 'Gray'), M, 0, 'Gray');
    end

    phix(n+1) = phix(n) - mu*imag(yx_hat(n)*conj(ephix));
    phiy(n+1) = phiy(n) - mu*imag(yy_hat(n)*conj(ephiy));
end

Y = [yx_hat; yy_hat];
phi = [phix; phiy];
