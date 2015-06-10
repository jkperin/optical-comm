% Calculate eigenvalues D and eigenfunctions Phi of the KL series expansion
% in the frequency domain.

% KWPhi = PhiD
% where K(vn, vm) = Ho(vn)He(vn-vm)Ho(vm), and W is a diagonal matrix of
% the weights of the Gauss-Legendre quadrature rule.

% Columns of Phi are orthornormal.

function [D, Phi, Fmax, nu] = klse_freq(rx, sim)

% Calculate Legendre-Gauss nodes and weights
H = rx.optfilt.H(sim.f/sim.fs);
H = cumtrapz(sim.f/sim.fs, abs(H).^2)/trapz(sim.f/sim.fs, abs(H).^2);
Fmax = sim.f(find(H >= 1-1e-4, 1, 'first'))/sim.fs;

% Legendre-Gauss quadrature
[nu,w] = lgwt(sim.Me, -Fmax, Fmax);
nu = nu(end:-1:1); % put in the right order
w = w(end:-1:1);

Ho = rx.optfilt.H; % annonymous functions of frequency response H(f/fs)
He = rx.elefilt.H;
Kf = @(f1, f2) Ho(f1).*He(f1-f2).*conj(Ho(f2));

% % Check accuracy
% ff = linspace(-Fmax, Fmax, 2^14);
% log(trapz(ff, abs(Ho(ff)).^2)) - log(sum(w.*abs(Ho(nu)).^2))
% log(trapz(ff, abs(He(ff)).^2)) - log(sum(w.*abs(He(nu).^2)))


K = zeros(sim.Me*[1 1]);
for k = 1:sim.Me
    K(:, k) = Kf(nu(k)*ones(sim.Me, 1), nu);
end

% Calculate eigenvalues
W = diag(w);
A = sqrt(W)*K*sqrt(W);
[B, D] = eig(A);

% D = real(diag(D)); % eigenvalues
D = real(diag(D));

% Rearrenge eigenvalues and eigenvectors in descending order
[D, ix] = sort(abs(D), 'descend');
B = B(:, ix);

Phi = diag(sqrt(1./w))*B; % eigenfunctions (not necessarily normalized)

% Normalize eigenfunctions and eigenvalues
for k = 1:sim.Me
    a = sum(w.*abs(Phi(:, k)).^2);
    
    Phi(:, k) = 1/sqrt(a)*Phi(:, k);
    
    D(k) = sqrt(a)*D(k);
    
    if sim.verbose && k <= 5
        figure(101)
        subplot(211), hold on, box on
        plot(nu, abs(Phi(:, k)).^2)
        xlabel('f/f_s', 'FontSize', 12)
        ylabel('|\phi_n(f/f_s)|^2', 'FontSize', 12)
        
        subplot(212), hold on, box on
        plot(nu, unwrap(angle(Phi(:, k))))
        xlabel('f/f_s', 'FontSize', 12)
        ylabel('phase', 'FontSize', 12)
    end
end



