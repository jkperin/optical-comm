function [Y, phi] = phase_estimation(X, eq, M, verbose)  
%% Adaptively calculates and correct for phase offset
Ntrain = eq.Ntrain; % Number of traning symbols
mu = eq.mu; % adaptation rate
Ytrain = eq.trainSeq; % training sequence 
Nsymb = length(X);

% Initialize filters coefficients
phix = zeros(1, Nsymb+1);
phiy = zeros(1, Nsymb+1);
ephix = zeros(1, Nsymb);
ephiy = zeros(1, Nsymb);

yx_hat = zeros(1, Nsymb); % output in x pol stream
yy_hat = zeros(1, Nsymb); % output in y pol stream
for n = 1:Nsymb % n runs over symbols
    % Phase tracker
    yx_hat(n) = X(1, n)*exp(-1j*phix(n));
    yy_hat(n) = X(2, n)*exp(-1j*phiy(n));

    if n < Ntrain % Training based on training sequence
        ephix(n) = yx_hat(n) - Ytrain(1, n);
        ephiy(n) = yy_hat(n) - Ytrain(2, n);
    else % based on deteced symbols
        ephix(n) = yx_hat(n) - qammod(qamdemod(yx_hat(n), M, 0, 'Gray'), M, 0, 'Gray');
        ephiy(n) = yy_hat(n) - qammod(qamdemod(yy_hat(n), M, 0, 'Gray'), M, 0, 'Gray');
    end

    phix(n+1) = phix(n) - mu*imag(yx_hat(n)*conj(ephix(n)));
    phiy(n+1) = phiy(n) - mu*imag(yy_hat(n)*conj(ephiy(n)));
end

Y = [yx_hat; yy_hat];
phi = [phix; phiy];

 % Phase tracker ----------------------------------------------------------
if exist('verbose', 'var') && verbose
    Ntrain_window = min(Ntrain, Nsymb);
    figure(406), clf
    subplot(211), hold on, box on
    plot(phi.')
    a = axis; 
    plot(Ntrain_window*[1 1], a(3:4), ':k')
    axis([1 Nsymb a(3:4)]);
    legend('X pol', 'Y pol', 'Training window')
    ylabel('Correction phase')
    xlabel('Symbol')
    title('Phase correction')

    subplot(212), hold on, box on
    plot(abs(ephix).^2)
    plot(abs(ephiy).^2)
    a = axis;
    plot(Ntrain_window*[1 1], a(3:4), ':k')
    axis([1 Nsymb a(3:4)]);
    legend('X pol', 'Y pol', 'Training window')
    ylabel('Mean square error')
    xlabel('Symbol')
    title('Mean square error')
end