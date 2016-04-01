function [Y, Wx, Wy, Wmix, MSE] = fse_cma2(eq, X)
%% Constant Modulus Algorithm
% Inputs:
% - eq : equalizer parameters {eq.ros = oversampling ratio, eq.Ntrain = 
% number of traning symbols, eq.mu = adaptation rate}
% - X : input samples in the two polarizations [2 x N] at rate ros x Rs
% Output:
% - Y : Equalized symbols at rate 1 x Rs
% - Wx, Wy : Filters coefficients 
% - MSE : mean square error
ros = eq.ros;
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
Wmix = zeros(2, Nfilters);  % crossover filters

% Initialize filters as if there were no ISI and there was no mixing between the two pols
Wx(ceil(Ntaps/2), :) = 1; 
Wy(ceil(Ntaps/2), :) = 1;

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
    zx = X(1, k + window);
    zy = X(2, k + window); 
  
    % Filter output
    yx = zx*Wx(:, filt);
    yy = zy*Wy(:, filt);
    
    % Filter output
    yx_hat(n) = yx + yy*Wmix(1, filt);
    yy_hat(n) = yx*Wmix(2, filt) + yy;

    % Calculate error
    ex(n) = 2 - abs(yx_hat(n))^2;
    ey(n) = 2 - abs(yy_hat(n))^2;

    % Update filters coefficients
    Wx(:, filt) = Wx(:, filt) + mu*ex(n)*yx_hat(n)*zx';
    Wy(:, filt) = Wy(:, filt) + mu*ey(n)*yy_hat(n)*zy';
    
    Wmix(1, filt) = Wmix(1, filt) + mu*ex(n)*yx_hat(n)*conj(yy);    
    Wmix(2, filt) = Wmix(2, filt) + mu*ey(n)*yy_hat(n)*conj(yx);    
    
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
