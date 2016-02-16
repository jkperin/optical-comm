function [Y, Wx, Wy, MSE] = adaptive_td_fse2(eq, X, Ytrain, M)
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
Ntrain = eq.Ntrain;
mu = eq.mu;

if mod(eq.Ntaps, 2) == 0 % always uses odd number of taps
    eq.Ntaps = eq.Ntaps + 1;
end
Ntaps = eq.Ntaps;

% Number of filters
[Nsamp, Nfilters]  = rat(ros);

% Initialize filters coefficients
Wx = zeros(Ntaps, Nfilters); % used for x polarization
Wy = zeros(Ntaps, Nfilters); % used for y polarization

% Initialize filters as if there were no ISI and there was no mixing between the two pols
Wx(ceil(Ntaps/2), :) = 1; 
Wy(ceil(Ntaps/2), :) = 1;
Wmix = eye(2); 

Nsymb = length(Ytrain);
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
    zx = X(1, k + window);
    zy = X(2, k + window); 

    % Filter output
    yx = zx*Wx(:, filt);
    yy = zy*Wy(:, filt);
       
    y = Wmix*[yx; yy];
    
    yx_hat(n) = y(1);
    yy_hat(n) = y(2);
    
    if n < Ntrain % Training based on training sequence
        if mod(n, 10) == 0
            ex(n) = yx - (Ytrain(1, n) - Wmix(1, 2)*yy)/Wmix(1, 1);
            ey(n) = yy - (Ytrain(2, n) - Wmix(2, 1)*yx)/Wmix(2, 2);
        else
            ex(n) = yx_hat(n) - Ytrain(1, n);
            ey(n) = yy_hat(n) - Ytrain(2, n);
        end
    else % based on deteced symbols
        ex(n) = yx_hat(n) - qammod(qamdemod(yx_hat(n), M), M);
        ey(n) = yy_hat(n) - qammod(qamdemod(yy_hat(n), M), M);
    end

    % Update filters coefficients
    if mod(n, 10) == 0
        Wx(:, filt) = Wx(:, filt) - 2*mu*ex(n)*zx';
        Wy(:, filt) = Wy(:, filt) - 2*mu*ey(n)*zy';
    else
        Wmix(1, :) = Wmix(1, :) - 2*mu*ex(n)*y.';
        Wmix(2, :) = Wmix(2, :) - 2*mu*ey(n)*y.';
    end
    
    if any(isnan(Wx))
        1;
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
1;