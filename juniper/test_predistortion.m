%% Test predistortion ideas
clear, clc

addpath f/
addpath ../f/
addpath ../mpam/

Nsymb = 2^12;
Mct = 10;
ros = 2;
N = Mct*Nsymb;
lamb = 1550e-9;

Fiber = fiber(20e3);

pulse_shape = select_pulse_shape('rect', ros);
mpam = PAM(4, 56e9, 'equally-spaced', pulse_shape);

fs = mpam.Rs*Mct;
[f, t] = freq_time(N, fs);
sim.Nsymb = Nsymb;
sim.Ndiscard = 31;
sim.f = f;
sim.fs = fs;
sim.Mct = Mct;

dataTX = randi([0 mpam.M-1], [1 Nsymb]);
% xt = 1/sqrt(2)*(xt + 1j*xt);
xd = sqrt(mpam.signal(dataTX));

xt = filter(ones(1, Mct/ros)*ros/Mct, 1, upsample(xd, Mct/ros));

% Perfect filtering
H = Fiber.Hdisp(f, lamb);
Eperfect = ifft(fft(xt).*ifftshift(conj(H).*exp(1j*2*pi*f/fs*(sim.Mct-1)/2)));

% optimized filter
p = design_predistortion_filter(31, ros, mpam, Fiber, lamb, true);

% F = design_filter('bessel', 5, 1/2*mpam.Rs/(fs/2));
% iEperfect = imag(Eperfect);
% iEperfect = real(ifft(fft(iEperfect).*ifftshift(F.H(f/fs))));
% Eperfect = real(Eperfect) + 1j*iEperfect;

hi = real(ifft(ifftshift((H))));
hq = imag(ifft(ifftshift((H))));
% % hq = fftshift(imag(ifft(ifftshift((H)))));
% E = ifft(fft(xt).*fft(ifftshift(hi) - 1j*ifftshift(hq)));

xdf = filter(p, 1, xd);
% xdf = xdf + 1j*circshift(xdf, [0 -1]);
% xdf = xd;
E = filter(ones(1, Mct/ros)*ros/Mct, 1, upsample(xdf, Mct/ros));
% E = real(ifft(fft(hi).*fft(E))) - 1j*real(ifft(fft(hi).*fft(E).*ifftshift(exp(-1j*pi/2*sign(f)))));
% E = real(ifft(fft(hq).*fft(E).*ifftshift(exp(1j*pi/2*sign(f))))) - 1j*real(ifft(fft(hq).*fft(E)));
% E = real(ifft(fft(hi).*fft(E))) - 1j*real(ifft(fft(hq).*fft(E)));

% figure(1082), hold on
% plot(fftshift(hi))
% plot(fftshift(hq))
% plot(real(fftshift(ifft(fft(hi).*ifftshift(exp(1j*pi/2*sign(f)))))))
% plot(imag(fftshift(ifft(fft(hi).*ifftshift(exp(1j*pi/2*sign(f)))))))

Eperfect = Eperfect/sqrt(mean(abs(Eperfect).^2));
E = E/sqrt(mean(abs(E).^2));

Eperfect = Fiber.linear_propagation(Eperfect, f, lamb);
E = Fiber.linear_propagation(E, f, lamb);

Pperfect = abs(Eperfect).^2;
P = abs(E).^2;

Rx.ADC.ros = 2;
Rx.ADC.fs = mpam.Rs*Rx.ADC.ros;
Rx.ADC.filt = design_filter('butter', 5, (Rx.ADC.fs/2)/(fs/2)); % Antialiasing filter
Rx.ADC.ENOB = Inf; % effective number of bits. Quantization is only included if sim.quantiz = true and ENOB ~= Inf
Rx.ADC.rclip = 0;
Rx.ADC.timeRefSignal = Pperfect;

[yk, ~, ytf] = adc(P, Rx.ADC, sim);

eq.type = 'Adaptive TD-LE';
eq.Ntaps = 15;
eq.Ntrain = Inf;
eq.ros = 2;
eq.mu = 1e-2;
eq.trainSeq = dataTX;
eq.Ndiscard = [1 1]*128;
[yd, eq] = equalize(eq, yk, [], mpam, sim, true);

figure(1), hold on, plot(yd, 'o')
axis([1 Nsymb -0.2 1.2])

% figure, hold on, box on
% plot(xt)
% plot(Pperfect)
% plot(P)
% 
% figure
% eyediagram(Pperfect, 2*Mct)

% figure, 
% eyediagram(P, 2*Mct)

