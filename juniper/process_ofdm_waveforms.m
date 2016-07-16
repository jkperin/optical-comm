%% Process scope data
clear, clc, close all

file = 'data/waveforms/wave-DMT-QPSK-59-ros=1_preemph.h5';
dacfile = 'data/waveforms/dmt_4qam_ros=1_waveform_59_preemph';
Dac = load(dacfile);

%
fs = 80e9; % Receiver sampling rate
ros = Dac.sim.ros.txDSP; % DAC oversampling rate
ofdm_ros = Dac.ofdm_ros; % OFDM oversampling rate: number of samples per OFDM symbol sample

%% Load channels
hinfo =  hdf5info(file);
ch1 = h5read(file, '/Waveforms/Channel 1/Channel 1Data');
ch2 = h5read(file, '/Waveforms/Channel 2/Channel 2Data');
ch3 = h5read(file, '/Waveforms/Channel 3/Channel 3Data');
ch4 = h5read(file, '/Waveforms/Channel 4/Channel 4Data');

ch1 = double(ch1).';
ch2 = double(ch2).';
ch3 = double(ch3).';
ch4 = double(ch4).';

[~, t] = freq_time(length(ch1), fs);

%% Orthornomalization
ch1 = ch1 - mean(ch1);
ch1 = ch1/sqrt(mean(abs(ch1).^2));
ch2 = ch2 - mean(ch2);
ch2 = ch2/sqrt(mean(abs(ch2).^2));
ch3 = ch3 - mean(ch3);
ch3 = ch3/sqrt(mean(abs(ch3).^2));
ch4 = ch4 - mean(ch4);
ch4 = ch4/sqrt(mean(abs(ch4).^2));

ch2 = ch2 - mean(ch1.*ch2)*ch1;
ch4 = ch4 - mean(ch3.*ch4)*ch4;

X = ch1 + 1j*ch2;
Y = ch3 + 1j*ch4;

%% Frequency offset compensaion
% [~, foff_ind] = max(abs(fftshift(fft(X))));
% foff = f(foff_ind)
% 
% X = X.*exp(-1j*2*pi*foff*t);
% Y = Y.*exp(-1j*2*pi*foff*t);
% 
% figure, plot(f, abs(fftshift(fft(X).^2)))
% figure, plot(t, real(X), t, imag(X))

%
Xrx = abs(X).^2; % Get intensity in X-pol

%% Antialiasing filter and reduce noise
f = freq_time(length(Xrx), fs);
Filt = design_filter('butter', 5, 1.12*(Dac.ofdm.fc(end))/(fs/2)); % half of the sampling rate
Xfilt = real(ifft(fft(Xrx).*ifftshift(Filt.H(f/fs))));

%% Resample
% Resample in order to have number of samples per symbol specified by ros
[p, q] = rat(ofdm_ros*Dac.ofdm.fs/fs);
Xrx = resample(Xfilt, p, q); % resample to match symbol rate at the DAC

Xrx = abs(Xrx).^2; % Get intensity in X-pol

xref = Dac.xk; % reference signal

% Align received signal
[c, lags] = xcorr(Xrx, xref);
cmax = max(c);
threshold = 0.9;
ind = find(c >= cmax*threshold & lags > 1);
pos = lags(ind);

figure, plot(lags, abs(c).^2)

Xrx = Xrx(pos(1):end);
Npatterns = floor(length(Xrx)/length(xref));
N = Npatterns*length(xref);
Xrx = Xrx(1:N);

%% OFDM symbol detection
ofdm = Dac.ofdm;

% Remove zeros
y = Xrx;
y = reshape(y, length(xref), Npatterns);
Npad = Dac.Npad;
y(1:floor(Npad/2), :) = [];
y(end-ceil(Npad/2)+1:end, :) = [];

% % Symbol-Trate sampling
y = y(1:ofdm_ros:end, :);

eq = Dac.Rx.AdEq;
eq.trainSeq = repmat(eq.trainSeq, 1, Npatterns);


% for k = 1:Npatterns
%     Xn = ofdm.detect(y(:, k), eq, true);
% end

[Xn, AGCn, W] = ofdm.detect(reshape(y, 1, []), eq, true);

figure, box on
plot(ofdm.fc/1e9, 10*log10(abs(AGCn.*W).^2), 'o')
title('Preemphasis necessary')
xlabel('Subcarrier frequency (GHz)')
ylabel('Preemphasis necessary (dB)')

ofdm.dataTX = repmat(ofdm.dataTX, 1, Npatterns);
ofdm.countBER([1024 1024], true)
