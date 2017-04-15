function [Y, W, MSE] = equalize_mimo4x2(X, eq, ModFormat, verbose)
%% Adaptive time-domain fractionally spaced linear equalizer adapted using LMS algorithm
% Inputs:
% - eq : equalizer parameters {eq.ros = oversampling ratio, eq.Ntrain = 
% number of traning symbols, eq.mu = adaptation rate}
% - X : input samples in the two polarizations [2 x N] at rate ros x Rs
% - Xtrain : symbols for training length(Xtrain) >= Ntrain
% - verbose (optional, default=false): whether to plot equalizer frequency
% response and convergence
% Output:
% - Y : Equalized symbols at rate 1 x Rs
% - Wx, Wy : Filters coefficients 
% - MSE: mean square error
M = ModFormat.M;  % qam order
Ytrain = eq.trainSeq;
ros = eq.ros;
Ntrain = eq.Ntrain;
mu = eq.mu;

if mod(eq.Ntaps, 2) == 0 % always uses odd number of taps
    eq.Ntaps = eq.Ntaps + 1;
end
Ntaps = eq.Ntaps;

% Initialize filters coefficients
W1 = zeros(4, Ntaps);
W1(:, ceil(Ntaps/2)) = 1; % Initialize filters as if there were no ISI
W2 = zeros(4, Ntaps);
W2(:, ceil(Ntaps/2)) = 1; % Initialize filters as if there were no ISI

Nsymb = length(Ytrain);
window = (1:Ntaps) - ceil(Ntaps/2); % {-(Ntaps-1)/2, ..., 0, ..., (Ntaps-1)/2}.
kstart = ceil(Ntaps/(2*ros))*ros+1; 
k = kstart; % runs over samples
y1_hat = zeros(1, Nsymb);
y2_hat = zeros(2, Nsymb);
e1 = zeros(2, Nsymb);
e2 = zeros(2, Nsymb);
for n = ((kstart-1)/ros+1):Nsymb-ceil(Ntaps/2) % n runs over symbols
    % k indexes samples and n indexes symbols
    % Build input
    z = X(:, k + window);

    % Filter output;
    y1_hat(n) = 0;
    y1_hat(n) = 0;
    for i = 1:4
        y1_hat(n) = y1_hat(n) + W1(i, :)*z(i, :).';
        y2_hat(n) = y2_hat(n) + W2(i, :)*z(i, :).';
    end

    % Calculate error
    if n < Ntrain % Training based on training sequence
        e1(n) = y1_hat(n) - Ytrain(1, n);
        e2(n) = y2_hat(n) - Ytrain(2, n);
    else % based on deteced symbols
        e1(n) = y1_hat(n) - pammod(pamdemod(y1_hat(n), M, 0, 'Gray'), M, 0, 'Gray');
        e2(n) = y2_hat(n) - pammod(pamdemod(y2_hat(n), M, 0, 'Gray'), M, 0, 'Gray');
    end

    % Update filters coefficients
    for i = 1:4
        W1(i, :) = W1(i, :) - 2*mu*e1(n)*conj(z(i, :));
        W2(i, :) = W1(i, :) - 2*mu*e2(n)*conj(z(i, :));
    end
    
    k = k + 1;
end

W = {W1, W2};

% Build outputs
Y = [y1_hat; y2_hat];
MSE = [abs(e1).^2; abs(e2).^2];


% Equalizer frequency response and convergence-----------------------------
if exist('verbose', 'var') && verbose
    figure(102)
    subplot(221), hold on, box on
    plot(1:length(MSE), MSE(1, :))
    a = axis;
    plot([1 1]*eq.Ntrain, [a(3) a(4)], '--k')
    axis([1 length(MSE) a(3:4)])
    xlabel('Iteration'); ylabel('MSE')
    
    subplot(222), hold on, box on
    plot(1:length(MSE), MSE(2, :))
    a = axis;
    plot([1 1]*eq.Ntrain, [a(3) a(4)], '--k')
    axis([1 length(MSE) a(3:4)])
    xlabel('Iteration'); ylabel('MSE')
    
    subplot(223), hold on, box on
    for filt = 1:4
        [hx, w] = freqz(W1(filt, :), 1);
        hy = freqz(W2(filt, :), 1, w);
        hline(1) = plot(w/(2*pi), abs(hx).^2, '-');
        hline(2) = plot(w/(2*pi), abs(hy).^2, '--', 'Color', get(hline(1), 'Color'));
    end
    xlabel('Frequency')
    ylabel('Magnitude')
    legend(hline, {'X pol', 'Y pol'})
    title(sprintf('%s, %d taps', eq.type, eq.Ntaps))
    
    subplot(224), hold on, box on
    for filt = 1:4
        [hx, w] = freqz(W1(filt, :), 1);
        hy = freqz(W2(filt, :), 1, w);
        hline(1) = plot(w/(2*pi),  unwrap(angle(hx)), '-');
        hline(2) = plot(w/(2*pi),  unwrap(angle(hy)), '--', 'Color', get(hline(1), 'Color'));
    end    
    xlabel('Frequency')
    ylabel('Phase')
    legend('X pol', 'Y pol')
    title(sprintf('%s, %d taps', eq.type, eq.Ntaps))
end
