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
Fmax = sim.f(find(H >= 1-1e-5, 1, 'first'))/sim.fs;

if Fmax < 0.25
    warning(sprintf('klse_fourier: Fmax = %.2f. Oversampling ratio must be increased so that Fmax > 0.25, otherwise A matrix in the KLSE will not be correctly calculated', Fmax));
end

% Legendre-Gauss quadrature
[nu,w] = lgwt(sim.Me, -Fmax, Fmax);
nu = nu(end:-1:1); % put in the right order
w = w(end:-1:1);

% Annonymous functions of frequency response H(f/fs)
He = rx.elefilt.H;
Ho = rx.optfilt.H(nu); 

% !! Since the frequency response functions are in discrete time, Fmax <= 0.25,
% otherwise He(f) might be calculated up to [-1, 1]. A warning is issued if
% Fmax > 0.25. This was also changed in klse_freq.
K = zeros(sim.Me);
for k = 1:sim.Me
    K(:, k) = Ho(k)*He(nu(k)-nu).*conj(Ho);
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

%
disp('klse_fourier: ratio between first and last eigenvalues in dB')
DU = 20*log10(abs(D(1)/D(end)))

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



