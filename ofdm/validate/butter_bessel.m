clear, clc, close all

n = 5;
fc = 0.8; % fc/(fs/2), fsct = fsMct = 1
Mct = 8;

wc = pi*fc/Mct;
w = linspace(0,pi,1e3);

% Butter
[bbutter, abutter] = butter(n, wc, 's');

% Bessel
wo = 1.65*wc;
[bbes, abes] = besself(n, wo);

hbutter = freqs(bbutter, abutter, w);
hbes = freqs(bbes, abes, w);

plot(w/(2*pi), abs(hbutter).^2, w/(2*pi), abs(hbes).^2)


% Convert to DT
% Butter
wcw = 2*tan(pi*fc/(2*Mct));
[bbutter, abutter] = butter(n, wcw, 's');
sys = tf(bbutter, abutter);
sysz = c2d(sys, 1, 'tustin'); % convert into DT by bilinear transformation
[bbutter, abutter] = tfdata(sysz);
bbutter = cell2mat(bbutter);
abutter = cell2mat(abutter);

% Bessel
wow = 2*tan(pi*1.65*fc/(2*Mct));
[bbes, abes] = besself(n, wow);
sys = tf(bbes, abes);
sysz = c2d(sys, 1, 'tustin'); % convert into DT by bilinear transformation
[bbes, abes] = tfdata(sysz);
bbes = cell2mat(bbes);
abes = cell2mat(abes);

[bbutter2, abutter2] = butter(n, fc/Mct);
hbutter2 = freqz(bbutter2, abutter2, w/(2*pi), 1);

hbutter = freqz(bbutter, abutter, w/(2*pi), 1);
hbes = freqz(bbes, abes, w/(2*pi), 1);

figure
plot(w/(2*pi), abs(hbutter).^2, w/(2*pi), abs(hbes).^2, w/(2*pi), abs(hbutter2).^2, '--r')


% wc = 2*tan(pi*rx.filter_cutoff/(2*sim.Mct)); % prewarped cutoff frequency
% wo = wc*sqrt(2^(1/rx.filter_order)-1);
% [bfilt, afilt] = besself(rx.filter_order, w0); % design prototype in CT with frequency prewarped
% sys = tf(bfilt, afilt);
% sysz = c2d(sys, 1, 'tustin'); % convert into DT by bilinear transformation
% [badc, aadc] = tfdata(sysz);
% badc = cell2mat(badc);
% aadc = cell2mat(aadc);
% 
% rx.gadc = impz(badc, aadc, rx.imp_length).';
% rx.gadc = rx.gadc/sum(rx.gadc);
% rx.gadc_delay = grpdelay(badc, aadc, 1);