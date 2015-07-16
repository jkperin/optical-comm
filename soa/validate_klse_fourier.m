%% Calulate MGF for the system with SOA
clear, clc, close all

addpath ../mpam/
addpath ../f
addpath f

% Simulation parameters
sim.M = 2; % Ratio of optical filter BW and electric filter BW
sim.Me = 16; 
sim.Mct = 17; % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
sim.Nsymb = 256;
sim.N = sim.Nsymb*sim.Mct;
sim.fs = mpam.Rs*sim.Mct;

sim.verbose = true;

% M-PAM
% M, Rb, leve_spacing, pshape
mpam = PAM(8, 100e9, 'equally-spaced', @(n) double(n >= 0 & n < sim.Mct));

% Time and Frequency
dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

sim.t = t;
sim.f = f;

f = f/sim.fs;

% Electric Lowpass Filter
rx.R = 1;
rx.elefilt = design_filter('bessel', 5, 1.25*mpam.Rs/(sim.fs/2));
He = rx.elefilt.H;

% Optical Bandpass Filter
rx.optfilt = design_filter('butter', 5, sim.M*rx.elefilt.fcnorm);
Ho = rx.optfilt.H;

% Generate optical signal
Nsymb = sim.Nsymb;
Nzero = 2*sim.Mct;
dataTX = randi([0 mpam.M-1], [Nsymb 1]);
Pd = pammod(dataTX, mpam.M, 0, 'gray');
pt = ones(1, sim.Mct);
Pt = reshape(kron(Pd, pt).', sim.N, 1);
Pt = Pt + (mpam.M - 1);
x = sqrt(Pt);
x([1:Nzero end-Nzero+1:end]) = 0;
X = fftshift(fft(x));

% Noise
% w = sqrt(soa.Seq(soa.G)*sim.fs/2)*(randn(sim.N, 1) + 1j*randn(sim.N, 1));

% et = x + w;
et = x; % noiseless
Ef = fftshift(fft(et));

eo = ifft(fft(et).*ifftshift(Ho(f)));

yt = real(ifft(fft(abs(eo).^2).*ifftshift(He(f))));

% figure
% subplot(211)
% plot(f, abs(He(f)).^2, f, abs(Ho(f)).^2)
% legend('electric', 'optical')
% 
% subplot(212)
% plot(f, unwrap(angle(He(f))), f, unwrap(angle(Ho(f))))
% legend('electric', 'optical')

% KL series expansion in the frequency domain
[U, D, Fmax] = klse_fourier(rx, sim, sim.N);

% Fourier series coefficients for periodic extension of x
xn = fftshift(fft(x))/length(x);
xn = xn(abs(f) <= Fmax);

fm = f(abs(f) <= Fmax);

% Calculate output
xk = zeros(length(fm), sim.N);
yk = zeros(length(x), 1); % output signal (noiseless)
for k = 1:length(x)
    xk(:, k) = (xn.*exp(1j*2*pi*fm*(k-1)));
    yk(k) = real(xk(:, k)'*(U*diag(D)*U')*xk(:, k));
end

% Used to calculate non-centrality parameter of Chi-Squared distributions
% ck(i, j) = ith chi-square distribution at jth time sample. |ck(i, j)|^2
% is the non-centrality parameter of that chi-square distribution
ck = U'*xk; 

figure, hold on
plot(yt)
plot(Pt)
plot(yk, '--k')
legend('numerical', 'Transmitted power', 'KL-SE Fourier')

% Discard zeros and downsample
ix = Nzero+sim.Mct/2:sim.Mct:sim.N-Nzero;
yd = yk(ix);
ck = ck(:, ix);
x = x(ix);

% plot
plot(ix, yd, 'o')

% SOA
% soa = soa(10^(10/10), 9, 1310e-9, 20);
% varASE = 1e5*(soa.N0*sim.fs/2)/sim.N;
% varTher = 0;
% 
% px = pdf_saddlepoint_approx(linspace(0, 2*max(Pt)), D, ck, varASE, varTher, true);


% Sampling
% yd = yt(sim.Mct/2:sim.Mct:end);
% Pd = Pt(sim.Mct/2:sim.Mct:end);
% yyd = yy(sim.Mct/2:sim.Mct:end);

% en = zeros(sim.N/sim.Mct, sim.Me);
% yyd= 0;
% fd = f(1:sim.Mct:end);
% for k = 1:sim.Me
%     phi = Phi(:, k);
%     phi = spline(nu, phi, fd/sim.fs);
%     en(:, k) = ifft(ifftshift(Ef(1:sim.Mct:end).*conj(phi)));
%    
%     yyd = yyd + D(k).*abs(en(:, k)).^2;
% end
% 
% figure, hold on
% stem(yd)
% stem(Pd)
% stem(yyd, 'x')
% legend('numerical', 'Transmitted power', 'KL-SE')


    