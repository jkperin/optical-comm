%% Analysi of Stokes vector receiver
clear, clc, close all

addpath ../f/

Nrealizations = 1e4;

P = zeros(Nrealizations, 4);
for k = 1:Nrealizations
    % Generate random unitary matrix
    phi = rand(1, 3)*2*pi;
    U1 = [exp(-1j*phi(1)/2), 0; 0 exp(1j*phi(1)/2)];
    U2 = [cos(phi(2)/2) -1j*sin(phi(2)/2); -1j*sin(phi(2)/2) cos(phi(2)/2)];
    U3 = [cos(phi(3)/2) -sin(phi(3)/2); sin(phi(3)/2) cos(phi(3)/2)];

    U = U1*U2*U3; % = [a -b; b^* a^*], for a and b complex

    a = U(1, 1);
    b = -U(1, 2);
%     
%     U*U'
%     
%     abs(a)^2 + abs(b)^2

    % Define fiber transfer matrix in Muller space
    % As defined in [1] M. Chagnon et al,"Digital signal processing for
    % dual-polarization intensity and interpolarization phase modulation 
    % formats using stokes detection", J. Lightw. Technol. 34, 188–195. 2016.
    M = [abs(a)^2, abs(b)^2, -2*real(a*b'), 2*imag(a*b');...
        abs(b)^2, abs(a)^2,  2*real(a*b'), -2*imag(a*b');...
        real(a*b), -real(a*b),  real(a^2) - real(b^2), -imag(a^2) - imag(b^2);...
        imag(a*b), -imag(a*b), imag(a^2) - imag(b^2), real(a^2) + real(b^2)];

    Mstar = inv(M);
    
    S(:, k) = svd(M);
        
    P(k, :) = diag(Mstar*Mstar.');
end

figure, hold on, box on
plot(P(:, 1), 'o')
plot(P(:, 2), 'x')
plot(P(:, 3), 's')
plot(P(:, 4), 'v')


figure, hold on, box on
boxplot(P)

% m = matlab2tikz(gca);
% m.write('stokes-noise-analysis.tex')

figure, 
subplot(221)
hist(P(:, 1), 50)
subplot(222)
hist(P(:, 2), 50)
subplot(223)
hist(P(:, 3), 50)
subplot(224)
hist(P(:, 4), 50)


