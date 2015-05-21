% KL series expansion in the frequency domain
function [xnd, D, yyd] = klse_freq(X, Ho, He, sim, Fmax)

f = sim.f;
t = sim.t;

% Calculate Legendre-Gauss nodes and weights
[nu,w] = lgwt(sim.Me, -Fmax, Fmax);

Kf = @(f1, f2) Ho(f1).*He(f1-f2).*Ho(f2);

K = zeros(sim.Me*[1 1]);
for k = 1:sim.Me
    K(:, k) = Kf(sim.fs*nu(k)*ones(sim.Me, 1), sim.fs*nu);
end

% Solve eigenvalue problem
W = diag(w);
A = sqrt(W)*K*sqrt(W);
[B, D] = eig(A);

D = real(diag(D)); % eigenvalues
Phi = sqrt(diag(1./w))*B; % orthonormal eigenfunctions

% Normalize eigenfunctions and eigenvalues
fs = f(1:sim.N/sim.Me:end)/sim.fs;
for k = 1:sim.Me
    a = trapz(fs, abs(Phi(:, k)).^2);
    
    Phi(:, k) = 1/sqrt(a)*Phi(:, k);
    
    D(k) = sqrt(a)*D(k);
end

en = zeros(sim.N, sim.Me);
xn = zeros(sim.N, sim.Me);
xnd = zeros(sim.Nsymb+2*sim.Nzero, sim.Me);
yy = 0;
yyd = 0;
for k = 1:sim.Me
    phi = Phi(:, k);
%     plot(nu, phi)

    % Interpolate before converting to time domain
    phi = spline(nu, phi, f/sim.fs);
    
%     if sim.verbose && k <= 5 
%         figure(100), hold on
%         plot(t, abs(phi).^2)
%         xlabel('Time (s)', 'FontSize', 12)
%         ylabel('|\phi_n(f)|^2', 'FontSize', 12)
%     end
    
    % Renormalize
    a = trapz(f/sim.fs, abs(phi).^2);
    phi = 1/sqrt(a)*phi;
    D(k) = sqrt(a)*D(k);
    
    % xn in 'continuous' time
    xn(:, k) = ifft(ifftshift(X.*conj(phi)));

%     en(:, k) = ifft(ifftshift(Ef.*conj(phi)));
%     yy = yy + D(k).*abs(en(:, k)).^2;
    
    % Xn sampled at decision instants (symbol rate)
    xnd(:, k) = xn(sim.Mct/2:sim.Mct:end, k);
    yyd = yyd + D(k).*abs(xnd(:, k)).^2;
end

% Remove zeros from begining and end of the sequence
nzero = [1:sim.Mct*sim.Nzero sim.N-sim.Mct*sim.Nzero+1:sim.N];
yyd([1:sim.Nzero sim.N/sim.Mct-sim.Nzero+1:sim.N/sim.Mct]) = [];
xnd([1:sim.Nzero sim.N/sim.Mct-sim.Nzero+1:sim.N/sim.Mct], :) = [];

