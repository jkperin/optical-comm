%% Calulate MGF for the system with SOA
clear, clc, close all

addpath f

% M-PAM
mpam.M = 4;
mpam.Rb = 100e9;
mpam.Rs = mpam.Rb/log2(mpam.M);
sim.Me = 64;

% Simulation parameters
sim.M = 1; % Ratio of optical filter BW and electric filter BW (must be integer)
sim.Mct = 8;
sim.N = 512;
sim.fs = mpam.Rs*sim.Mct;

dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

% Electric Lowpass Filter
EleFilt.type = 'gaussian';
EleFilt.order = 4;
EleFilt.BW = mpam.Rs; % single-sided BW
[be, ae] = design_filter(EleFilt.type, EleFilt.order, EleFilt.BW/(sim.fs/2), sim.Mct);
EleFilt.grpdelay = grpdelay(be, ae, 1);
He = @(f) freqz(be, ae, f, sim.fs).*exp(1j*2*pi*f/sim.fs*EleFilt.grpdelay);
% He = @(f) freqz(be, ae, f, sim.fs);

% Optical Bandpass Filter
OptFilt.type = 'gaussian';
OptFilt.order = 4;
OptFilt.BW = sim.M*EleFilt.BW; % single-sided BW
[bo, ao] = design_filter(OptFilt.type, OptFilt.order, OptFilt.BW/(sim.fs/2), sim.Mct);
OptFilt.grpdelay = grpdelay(bo, ao, 1);
Ho = @(f) freqz(bo, ao, f, sim.fs).*exp(1j*2*pi*f/sim.fs*OptFilt.grpdelay);
% Ho = @(f) freqz(bo, ao, f, sim.fs);
% 
tx.lamb = 1310e-9;
soa.Fn = 9; % dB
soa.Seq = @(G) (G - 1)*10^(soa.Fn/10)/2*(1.98644568e-25/tx.lamb); % one-sided PSD
soa.G = 10^(5/10);

% Generate optical signal
Nzero = 32;
Ptx = 1e-3;
Nsymb = sim.N/sim.Mct;
dataTX = randi([0 mpam.M-1], [Nsymb 1]);
Pd = pammod(dataTX, mpam.M, 0, 'gray');
pt = ones(1, sim.Mct);
Pt = reshape(kron(Pd, pt).', sim.N, 1);
Pt = Pt + (mpam.M - 1);
Pt = Pt/mean(Pt)*Ptx*soa.G;
Pt(1:Nzero) = 0;
Pt(end-Nzero+1:end) = 0;
x = sqrt(Pt);
X = fftshift(fft(x));

% Noise
w = 0*sqrt(soa.Seq(soa.G)*sim.fs/2)*(randn(sim.N, 1) + 1j*randn(sim.N, 1));

et = x + w;
Ef = fftshift(fft(et));

eo = ifft(fft(et).*ifftshift(Ho(f)));

yt = real(ifft(fft(abs(eo).^2).*ifftshift(He(f))));

figure
plot(f/sim.fs, abs(He(f)).^2, f/sim.fs, abs(Ho(f)).^2)
legend('electric', 'optical')

% KL series expansion in the frequency domain
Fmax = 0.5; %1.5*OptFilt.BW/sim.fs;

[nu,w] = lgwt(sim.Me, -Fmax, Fmax);

Kf = @(f1, f2) Ho(f1).*He(f1-f2).*Ho(f2);

K = zeros(length(nu)*[1 1]);
for k = 1:length(nu)
    K(:, k) = Kf(sim.fs*nu(k)*ones(length(nu), 1), sim.fs*nu);
end

W = diag(w);

A = sqrt(W)*K*sqrt(W);

[B, D] = eig(A);
D = real(diag(D));

Phi = sqrt(diag(1./w))*B;

en = zeros(sim.N, sim.Me);
yy = 0;
for k = 1:sim.Me
    phi = Phi(:, k);
%     plot(nu, phi)
    phi = spline(nu, phi, f/sim.fs);
    en(:, k) = ifft(ifftshift(Ef.*conj(phi)));
   
    yy = yy + D(k).*abs(en(:, k)).^2;
end

yt = (yt-mean(yt))/std(yt);
Pt = (Pt-mean(Pt))/std(Pt);
yy = (yy - mean(yy))/std(yy);

figure, hold on
plot(yt)
plot(Pt)
plot(yy, '--k')

    