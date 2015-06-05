%% Calulate MGF for the system with SOA
clear, clc, close all

addpath f

% M-PAM
mpam.M = 4;
mpam.Rb = 100e9;
mpam.Rs = mpam.Rb/log2(mpam.M);
sim.Me = 64; % Number of eigenvalues
sim.verbose = true;

% Simulation parameters
sim.M = 2; % Ratio of optical filter BW and electric filter BW (must be integer)
sim.Mct = 8;
sim.N = 512;
sim.fs = mpam.Rs*sim.Mct;

dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';
f = f/sim.fs;

sim.t = t;
sim.f = f;

% Electric Lowpass Filter
rx.elefilt = design_filter('bessel', 5, 1.25*mpam.Rs/(sim.fs/2));
He = rx.elefilt.H;

% Optical Bandpass Filter
rx.optfilt = design_filter('bessel', 5, sim.M*rx.elefilt.fcnorm);
Ho = rx.optfilt.H;

% Generate optical signal
Nsymb = sim.N/sim.Mct;
dataTX = randi([0 mpam.M-1], [Nsymb 1]);
Pd = pammod(dataTX, mpam.M, 0, 'gray');
pt = ones(1, sim.Mct);
Pt = reshape(kron(Pd, pt).', sim.N, 1);
Pt = Pt + (mpam.M - 1);
x = sqrt(Pt);
X = fftshift(fft(x));

% Noise
% w = sqrt(soa.Seq(soa.G)*sim.fs/2)*(randn(sim.N, 1) + 1j*randn(sim.N, 1));

% et = x + w;
et = x; % noiseless
Ef = fftshift(fft(et));

eo = ifft(fft(et).*ifftshift(Ho(f)));

yt = real(ifft(fft(abs(eo).^2).*ifftshift(He(f))));

figure
subplot(211)
plot(f, abs(He(f)).^2, f, abs(Ho(f)).^2)
legend('electric', 'optical')

subplot(212)
plot(f, unwrap(angle(He(f))), f, unwrap(angle(Ho(f))))
legend('electric', 'optical')

% KL series expansion in the frequency domain
[D, Phi, Fmax, nu] = klse_freq(rx, sim);

D

en = zeros(sim.N, sim.Me);
yy = 0;
for k = 1:sim.Me
    phi = Phi(:, k);
%     plot(nu, phi)
    phi = spline(nu, phi, f);
    
    % Rescale after interpolation
    a = trapz(f, abs(phi).^2);
    phi = 1/sqrt(a)*phi;
    
    if k <= 3
        figure(100), hold on, box on
        plot(f, abs(phi).^2)
        xlabel('f/f_s', 'FontSize', 12)
        ylabel('|\phi_n(f/f_s)|^2', 'FontSize', 12)
    end
    legend('n = 1', 'n = 2', 'n = 3', 'n = 4', 'n = 5')

       
    en(:, k) = ifft(ifftshift(Ef.*conj(phi)));
   
    yy = yy + D(k).*abs(en(:, k)).^2;
end

figure, hold on
plot(yt)
plot(Pt)
plot(yy, '--k')
legend('numerical', 'Transmitted power', 'KL-SE')

% Sampling
yd = yt(sim.Mct/2:sim.Mct:end);
Pd = Pt(sim.Mct/2:sim.Mct:end);
yyd = yy(sim.Mct/2:sim.Mct:end);

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

figure, hold on
stem(yd)
stem(Pd)
stem(yyd, 'x')
legend('numerical', 'Transmitted power', 'KL-SE')


    