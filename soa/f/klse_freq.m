% Calculate eigenvalues D and eigenfunctions Phi of the KL series expansion
% in the frequency domain.

% KWPhi = PhiD
% where K(vn, vm) = Ho(vn)He(vn-vm)Ho(vm), and W is a diagonal matrix of
% the weights of the Gauss-Legendre quadrature rule.

% Columns of Phi are orthornormal.

function [D, Phi, Fmax, nu] = klse_freq(rx, sim)

% Calculate Legendre-Gauss nodes and weights
Fmax = min(1.5*rx.optfilt.fcnorm/2, 0.5) % Maximum frequency (same as used in ber_soa_klse_freq.m)

% Legendre-Gauss quadrature
[nu,w] = lgwt(sim.Me, -Fmax, Fmax);

sum(w.*abs(rx.optfilt.H(nu)).^2)/rx.optfilt.noisebw(1)

Ho = rx.optfilt.H; % annonymous functions of frequency response H(f/fs)
He = rx.elefilt.H;
Kf = @(f1, f2) Ho(f1).*He(f1-f2).*Ho(f2);

K = zeros(sim.Me*[1 1]);
for k = 1:sim.Me
    K(:, k) = Kf(nu(k)*ones(sim.Me, 1), nu);
end

% Calculate eigenvalues
W = diag(w);
A = sqrt(W)*K*sqrt(W);
[B, D] = eig(A);

D = real(diag(D)); % eigenvalues

% Rearrenge eigenvalues and eigenvectors in descending order
% [D, ix] = sort(D, 'descend');
% B = B(:, ix);

Phi = diag(sqrt(1./w))*B; % eigenfunctions (not necessarily normalized)

% Normalize eigenfunctions and eigenvalues
for k = 1:sim.Me
    a = sum(w.*abs(Phi(:, k)).^2);
    
    Phi(:, k) = 1/sqrt(a)*Phi(:, k);
    
    D(k) = sqrt(a)*D(k);
    
    if sim.verbose && k <= 5
        figure(100), hold on, box on
        plot(nu, abs(Phi(:, k)).^2)
        xlabel('f/f_s', 'FontSize', 12)
        ylabel('|\phi_n(f/f_s)|^2', 'FontSize', 12)
    end
end



