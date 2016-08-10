%% Validate OFDM waveform generation and processing
clear, clc, close all

dacfile = 'data/waveforms/ofdm4_Rb=42Gbps_19';

fs = 80e9; % DSO sampling rate
Dac = load(dacfile);

xq = Dac.xq; % quantized signal at the DAC sampling rate

[p, q] = rat(fs/Dac.Tx.DAC.fs);

x = resample(xq, p, q); % signal at the DSO sampling frequency

x = sin(pi/2*x/255); % modulator
x = repmat(x, 1, 100); % repeat
x = circshift(x, [0 2364]);
x = x.';

WaveForms = zeros(length(x), 4);
WaveForms(:, 1) = x;
WaveForms(:, 2) = 0*1e-3*randn(size(x));
WaveForms(:, 3) = 0*0.7*x + 0*1e-3*randn(size(x));
WaveForms(:, 4) = 0*1e-3*randn(size(x));

ber_count = process_ofdm_waveforms(WaveForms, dacfile);
