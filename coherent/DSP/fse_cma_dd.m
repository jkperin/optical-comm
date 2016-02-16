function [Y, Wx, Wy, MSE] = fse_cma_dd(eq, X, M)
%% Adaptive time-domain fractionally spaced linear equalizer adapted using LMS algorithm
% Inputs:
% - eq : equalizer parameters {eq.ros = oversampling ratio, eq.Ntrain = 
% number of traning symbols, eq.mu = adaptation rate}
% - X : input samples in the two polarizations [2 x N] at rate ros x Rs
% - Xtrain : symbols for training length(Xtrain) >= Ntrain
% - M : qam order
% - MSE: mean square error
% Output:
% - Y : Equalized symbols at rate 1 x Rs
% - Wx, Wy : Filters coefficients 
ros = eq.ros;
mu = eq.mu;
Ntrain = eq.Ntrain;

if mod(eq.Ntaps, 2) == 0 % always uses odd number of taps
    eq.Ntaps = eq.Ntaps + 1;
end
Ntaps = eq.Ntaps;

% Number of filters
[Nsamp, Nfilters]  = rat(ros);

% Initialize filters coefficients
Wx = zeros(2*Ntaps, Nfilters); % used for x polarization
Wy = zeros(2*Ntaps, Nfilters); % used for y polarization

% Initialize filters as if there were no ISI and there was no mixing between the two pols
Wx(ceil(Ntaps/2), :) = 1; 
Wy(ceil(Ntaps/2)+Ntaps, :) = 1;

Nsymb = floor(length(X)/ros);
ex = zeros(1, Nsymb);
ey = zeros(1, Nsymb);
yx_hat = zeros(1, Nsymb); % output in x pol stream
yy_hat = zeros(1, Nsymb); % output in y pol stream
filt = 1;
window = (1:Ntaps) - ceil(Ntaps/2); % {-(Ntaps-1)/2, ..., 0, ..., (Ntaps-1)/2}.
kstart = ceil(Ntaps/(2*Nsamp))*Nsamp+1; 
k = kstart; % runs over samples
for n = ((kstart-1)/ros+1):Nsymb-ceil(Ntaps/2) % n runs over symbols
    % k indexes samples and n indexes symbols
    z = [X(1, k + window), X(2, k + window)]; % input of filter [x pol, y pol]

    % Filter output
    yx_hat(n) = z*Wx(:, filt);
    yy_hat(n) = z*Wy(:, filt);

    % Calculate error
    if n < Ntrain % blind equalization
        ex(n) = 2 - abs(yx_hat(n))^2;
        ey(n) = 2 - abs(yy_hat(n))^2;
        
        % Update filters coefficients
        Wx(:, filt) = Wx(:, filt) + mu*ex(n)*yx_hat(n)*z';
        Wy(:, filt) = Wy(:, filt) + mu*ey(n)*yy_hat(n)*z';
    else % decision-directed error
        ex(n) = yx_hat(n) - qammod(qamdemod(yx_hat(n), M, 0, 'Gray'), M, 0, 'Gray');
        ey(n) = yy_hat(n) - qammod(qamdemod(yx_hat(n), M, 0, 'Gray'), M, 0, 'Gray');
        
        % Update filters coefficients
        Wx(:, filt) = Wx(:, filt) - 2*mu*ex(n)*z';
        Wy(:, filt) = Wy(:, filt) - 2*mu*ey(n)*z';
    end
    
    % Increment filter index
    filt = mod(filt, Nfilters) + 1;

    % ensures that the center of the window of samples remains 
    % close to the nth symbol
    if abs((k+1)/ros - n-1) > 0.5
        k = k + 2;
    else
        k = k + 1;
    end
end

% Build outputs
Y = [yx_hat; yy_hat];
MSE = [abs(ex).^2; abs(ey).^2];

figure(100), subplot(211)
plot(MSE(1, :))
subplot(212)
plot(MSE(2, :))
xlabel('Iteration')
ylabel('MSE')
