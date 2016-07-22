function ber_count = process_ofdm_waveforms(file_or_data, dacfile)
%% Process PAM waveforms saved with DSO

addpath f/
addpath ../f/
addpath ../mpam/

% file = 'data/waveforms/wave-DMT-QPSK-59-ros=1_preemph.h5';
% dacfile = 'data/waveforms/ofdm16_Rb=56Gbps_51';
Dac = load(dacfile);
ofdm = Dac.ofdm;

%
fs = 80e9; % Receiver sampling rate

%% Load channels
if ischar(file_or_data)
    file = file_or_data;
    ch1 = h5read(file, '/Waveforms/Channel 1/Channel 1Data');
    ch2 = h5read(file, '/Waveforms/Channel 2/Channel 2Data');
    ch3 = h5read(file, '/Waveforms/Channel 3/Channel 3Data');
    ch4 = h5read(file, '/Waveforms/Channel 4/Channel 4Data');
elseif strcmpi(class(file_or_data), 'double')
    data = file_or_data;
    ch1 = data(:, 1);
    ch2 = data(:, 2);
    ch3 = data(:, 3);
    ch4 = data(:, 4);
else
    error('process_pam_waveforms: First parameter must be either a .h5 file path or a matrix with the waveforms in each column')
end

ch1 = double(ch1).';
ch2 = double(ch2).';
ch3 = double(ch3).';
ch4 = double(ch4).';

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

%% Pre-filtering
% ch1 = real(ifft(fft(ch1).*ifftshift(Filt.H(f/fs))));
% ch2 = real(ifft(fft(ch2).*ifftshift(Filt.H(f/fs))));
% ch3 = real(ifft(fft(ch3).*ifftshift(Filt.H(f/fs))));
% ch4 = real(ifft(fft(ch4).*ifftshift(Filt.H(f/fs))));

% X = ch1 + 1j*ch2;
% Y = ch3 + 1j*ch4;

%% Convert to intensity
Xrx = abs(X).^2 + abs(Y).^2; % Get intensity in X-pol

%% Antialiasing filter and reduce noise
f = freq_time(length(Xrx), fs);
Filt = design_filter('butter', 5, 1.12*(Dac.ofdm.fc(end))/(fs/2)); % half of the sampling rate
Xfilt = real(ifft(fft(Xrx).*ifftshift(Filt.H(f/fs))));

Xfilt = Xfilt - mean(Xfilt);
Xfilt = Xfilt/sqrt(mean(abs(Xfilt).^2));

%% Resample
% Resample to DAC sampling rate
[p, q] = rat(Dac.Tx.DAC.fs/fs);
Xrx = resample(Xfilt, p, q); % resample to match symbol rate at the DAC

xref = Dac.xkDAC; % reference signal at OFDM sampling rate

% Align received signal
[c, lags] = xcorr(Xrx, xref);
cmax = max(c);
threshold = 0.95;
ind = find(c >= cmax*threshold & lags >= 1);
% select only highest correlation peaks

pos = lags(ind);

figure(201), clf, hold on
plot(lags, abs(c).^2)
plot(pos, abs(c(ind)).^2, 'xk', 'MarkerSize', 6);
xlabel('Lags')
ylabel('Crosscorrelation')
title('Aligning received signal to transmitted sequence')

% window = 0:length(xref)-1;
% 
% pstart = pos(1);
% Npatterns = 1;
% y = zeros(length(xref), ceil((length(Xrx)-pos(1))/length(xref)));
% while pstart + window(end) <= length(Xrx)
%     y(:, Npatterns) = Xrx(pstart + window).';
%     pstart = pos(find(pos >= pstart + window(end) + 1, 1, 'first'));
%     Npatterns = Npatterns + 1;
% end

Xrx = Xrx(pos(1):end);
Npatterns = floor(length(Xrx)/length(xref));
N = Npatterns*length(xref);
Xrx = Xrx(1:N);
% Remove zeros
Xrx = reshape(Xrx, [], Npatterns);
Xrx([1:floor(Dac.Npad/2) end-ceil(Dac.Npad/2)+1:end], :) = [];

% Downsample to OFDM sampling rate
y = reshape(Xrx, 1, []);
[p, q] = rat(ofdm.fs/Dac.Tx.DAC.fs);
y = resample(y, p, q);

%% OFDM symbol detection
try % compatibility reasons
    eq = Dac.AdEq;
catch err
    eq =  Dac.Rx.AdEq;
end
eq.mu = 1e-2;
eq.trainSeq = repmat(eq.trainSeq, 1, Npatterns);

[Xn, AGCn, W] = ofdm.detect(y, eq, true);

ofdm.dataTX = repmat(ofdm.dataTX, 1, Npatterns);
ber_count = ofdm.countBER([1000 20], true)
