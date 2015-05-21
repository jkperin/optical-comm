%% Verify how well discrete-time filter design approximates continuous-time filter design
clear, clc, close all

Mct = 8;            % oversampling ratio to emulate CT
fs = 100e9;         % sampling rate before oversampling
f = linspace(0, Mct*fs/2);
fc = 0.8;           % Normalized cut-off frequency before oversampling (fc/(fs/2))

% Butterworth
[bz,az] = butter(5, fc/Mct); % design in discrete time
hz = freqz(bz,az, f, fs*Mct);

[bs,as] = butter(5, 2*fs*Mct*tan(pi*fc/(2*Mct)), 's'); 
hs = freqs(bs,as, 2*pi*f);

figure
subplot(211)
plot(f/1e9, 20*log10(abs(hs)), f/1e9, 20*log10(abs(hz)))
xlabel('Frequency (GHz)')
ylabel('Magnitude |H(j\omega)|^2')
legend('CT', 'DT')
axis([0 fs/1e9 -80 5])
title(sprintf('5th-order Butterworth (Oversampling = %d)', Mct))
subplot(212)
plot(f/1e9, rad2deg(unwrap(angle(hs))), f/1e9, rad2deg(unwrap(angle(hz))))
xlabel('Frequency (GHz)')
ylabel('Phase (deg)')
legend('CT', 'DT')
axis([0 fs/1e9 -400 0])

clear bz az hz bs as hs

% Chebyshev type 1 (ripple on the passband)
ripple = 0.1; 
[bz,az] = cheby1(5, ripple, fc/Mct);
hz = freqz(bz,az, f, fs*Mct);
[bs,as] = cheby1(5, ripple, 2*pi*fc*fs/2, 's');
hs = freqs(bs, as, 2*pi*f);

figure
subplot(211)
plot(f/1e9, 20*log10(abs(hs)), f/1e9, 20*log10(abs(hz)))
xlabel('Frequency (GHz)')
ylabel('Magnitude |H(j\omega)|^2')
legend('CT', 'DT')
axis([0 fs/1e9 -80 5])
title(sprintf('5th-order Chebyshev 1 (%.1f dB of ripple on the passband, Oversampling = %d)', ripple, Mct))
subplot(212)
plot(f/1e9, rad2deg(unwrap(angle(hs))), f/1e9, rad2deg(unwrap(angle(hz))))
xlabel('Frequency (GHz)')
ylabel('Phase (deg)')
legend('CT', 'DT')
axis([0 fs/1e9 -400 0])

clear bz az hz bs as hs

% Bessel filter
wc2wo = 1.65;
wc = 2*pi*fc*fs/2;
wo = wc*wc2wo;
[bs, as] = besself(5, wo);
hs = freqs(bs, as, 2*pi*f);
wow = 2*fs*Mct*tan(pi*wc2wo*fc/(2*Mct)); % wo warped
[bs2, as2] = besself(5, wow); % CT with frequency prewarped
sys = tf(bs2, as2);
sysz = c2d(sys, 1/(fs*Mct), 'tustin'); % convert into DT by bilinear transformation
[bz, az] = tfdata(sysz);
bz = cell2mat(bz);
az = cell2mat(az);
hz = freqz(bz,az,f,fs*Mct);

figure
subplot(211)
plot(f/1e9, 20*log10(abs(hs)), f/1e9, 20*log10(abs(hz)))
xlabel('Frequency (GHz)')
ylabel('Magnitude |H(j\omega)|^2')
legend('CT', 'DT')
axis([0 fs/1e9 -80 5])
title(sprintf('5th-order Besself (Oversampling = %d)', Mct))
subplot(212)
plot(f/1e9, rad2deg(unwrap(angle(hs))), f/1e9, rad2deg(unwrap(angle(hz))))
xlabel('Frequency (GHz)')
ylabel('Phase (deg)')
legend('CT', 'DT')
axis([0 fs/1e9 -400 0])

clear bz az hz bs as hs
