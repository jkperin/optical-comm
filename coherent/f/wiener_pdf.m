%% Wiener process distribution
clear, clc, close all

varp = 1;

N = 2^14;
Ntrials = 256;
theta = zeros(Ntrials, N);
thetaq = theta;
for n = 1:Ntrials
    initial_phase    = pi*(2*rand(1)-1); % [-pi, pi]
    dtheta      = [0 sqrt(varp)*randn(1, N-1)];                                  % i.i.d.Gaussian random variables with zero mean and variance sigma_p^2
    theta(n, :) = initial_phase + cumsum(dtheta, 2);
    thetaq(n, :) = asin(sin(theta(n, :)));
end

figure, histfit(theta(:, 100))
figure, histfit(theta(:, end))
var(theta(:, 100))


figure, histfit(thetaq(:, 100))
figure, histfit(thetaq(:, end))
