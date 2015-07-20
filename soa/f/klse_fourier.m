%% KLSE Fourier Series Expansion 
% Calculate eigenvectors and eigenvalues for KLSE Fourier Series Expansion 
% Both noise and signal are expanded on the same basis.
% N = sequence length
% Hdisp (optional) = dispersion frequency response
function [U, D, Fmax] = klse_fourier(rx, sim, N, Hdisp)
%% Expansion using same basis for noise and signal. 
% % Input signal x is considered periodic with period L = length(x). Both
% % noise and signal are expanded over the same Fourier Series basis.
% % Problem: !! has to calculate eigenvalues of L x L matrix
%
if nargin < 4 % if fiber is not included
    Hdisp = @(f) 1;
end
    
% Bandwidth to account for great part of energy
H = rx.optfilt.H(sim.f/sim.fs).*Hdisp(sim.f);
H = cumtrapz(sim.f/sim.fs, abs(H).^2)/trapz(sim.f/sim.fs, abs(H).^2);
Fmax = sim.f(find(H >= 1-1e-4, 1, 'first'))/sim.fs;

if Fmax < 0.25
    warning(sprintf('klse_fourier: Fmax = %.2f. Oversampling ratio must be increased so that Fmax > 0.25, otherwise A matrix in the KLSE will not be correctly calculated', Fmax));
end

% Redefine frequency
df = 1/N;
ft = (-0.5:df:0.5-df).';
f = ft(abs(ft) <= Fmax);
L = length(f);

% Annonymous functions of frequency response H(f/fs)
He = rx.elefilt.H;
Ho = rx.optfilt.H(f).*Hdisp(f*sim.fs); 

% !! Since the frequency response functions are in discrete time, Fmax <= 0.25,
% otherwise He(f) might be calculated up to [-1, 1]. A warning is issued if
% Fmax > 0.25. This was also changed in klse_freq.
A = zeros(L);
for k = 1:L
    A(:, k) = Ho(k)*He(f(k)-f).*conj(Ho);
end

% Include photodiode responsivity
A = rx.R*A;

% Calcualte eigenvalues
if sim.Me >= L - 1
    [U, D] = eig(A);
else
    [U, D] = eigs(A, sim.Me);
end

% check
% norm(A - U*D*U')

D = real(diag(D));

%
DU = 20*log10(abs(D(1)/D(end)));
% if abs(DU) < 20
%     warning('klse_fourier: ratio between first and last eigenvalues in dB = %f\n', DU)
% end

%% Expansion using different bases for noise and signal (reduces number
% of eigenvalues)
% % Fourier series expansion (V.A) in "On the Error Probability Evaluation in Lightwave 
% % Systems With Optical Amplification"
% % Receiver memory T0
% mu = 2;
% T0 = mu*(1/rx.optfilt.noisebw(1) + 1/rx.elefilt.noisebw(1));
% 
% % Bandwidth of 99.99% of filter energy
% H = rx.optfilt.H(sim.f);
% H = cumtrapz(sim.f, abs(H).^2)/trapz(sim.f, abs(H).^2);
% B0 = sim.f(find(H >= 0.99999, 1, 'first'));
% 
% % Number of eigenvalues
% M = ceil(B0*T0);
% Me = 2*M + 1;
% 
% % Annonymous functions of frequency response H(f/fs)
% Ho = rx.optfilt.H; 
% He = rx.elefilt.H;
% Kf = @(f1, f2) Ho(f1).*He(f1-f2).*Ho(f2);
% 
% A = zeros(Me);
% nu = ((-M:M)/T0).';
% for k = 1:Me
%     A(:, k) = Kf(nu, ones(Me, 1)*(k - M - 1)/T0);
% end
% 
% % Calcualte eigenvalues
% [U, D] = eig(A);
% 
% % Rearrenge eigenvalues and eigenvectors in descending order
% D = diag(D);
% [D, ix] = sort(D, 'descend');
% U = U(:, ix);
% 
% xn = (fftshift(fft(x)).')/length(x);
% df = 1/length(x);
% f = -0.5:df:0.5-df;
% L = length(f);
% vk = zeros(Me, L); % each column corresponds to a vector vk
% for i = 1:Me
%     vk(i, :) = L*ifft(ifftshift(Kf(f, ones(1, L)*(i-M-1)/T0).*xn));
% end
% 
% % bk
% bk = U'*vk;
% 
% dk = zeros(L, 1);
% dt = 1/L;
% 
% % Equivalent yk (noiseless)
% xp = zeros(L);
% for l = 1:L
%     K = Kf(f, ones(1, L)*f(l));
%     xp(l, :) = L*ifft(ifftshift(xn.*K));
% end
% 
% for k = 1:L
%     dk(k) = real(sum(xp(:, k).*(xn').*exp(-1j*2*pi*((-L/2:L/2-1)')*(k-1)*dt)));
% end
% 
% % Equivalent to this
% % for k = 1:L
% %     for m = 1:L
% %         for l = 1:L
% %             K = Kf(f(m)*[1 1], f(l)*[1 1]);
% %             dk(k) = dk(k) + xn(m)*xn(l)'*K(1)*exp(1j*2*pi*(m-l)*(k-1)*dt);
% %         end
% %     end
% % end

%% Expansion using same basis for noise and signal (incomplete)
% Alternative Expansion (V.B) in "On the Error Probability Evaluation in Lightwave 
% Systems With Optical Amplification"

% % Receiver memory T0
% mu = 2;
% T0 = mu*(1/rx.optfilt.noisebw(1) + 1/rx.elefilt.noisebw(1));
% 
% % Bandwidth of 99.99% of filter energy
% H = rx.optfilt.H(sim.f);
% H = cumtrapz(sim.f, abs(H).^2)/trapz(sim.f, abs(H).^2);
% B0 = sim.f(find(H >= 0.99999, 1, 'first'));
% 
% % Number of eigenvalues
% M = 2*ceil(T0);
% Me = 2*M + 1;
% 
% f = -M:M/T0;
% 
% % Annonymous functions of frequency response H(f/fs)
% Ho = rx.optfilt.H; 
% He = rx.elefilt.H;
% Kf = @(f1, f2) Ho(f1).*He(f1-f2).*Ho(f2);
% 
% A = zeros(L);
% for k = 1:L
%     A(:, k) = Kf(ones(L, 1)*f(k), f);
% end
% 
% % Calcualte eigenvalues
% % [U, D] = eigs(A, 32);
% [U, D] = eig(A);
% 
% 
% W = ceil(T0/2);
% L = length(x);
% dk = zeros(1, L);
% xk = zeros(Me, L);
% dt = T0/Me;
% for k = W+1:length(x)-W
%     xt = interp1(k-W:k+W, x(k-W:k+W), k-T0/2:dt:k+T0/2-dt);
%     xn = (fftshift(fft(xt)).')/Me;
%     xk(:, k) = xn.*exp(1j*2*pi*(-M:M).'*k/T0);
%     dk(k) = real(xk(:, k)'*A*xk(:, k));
% end
% 
% bk = 0;

