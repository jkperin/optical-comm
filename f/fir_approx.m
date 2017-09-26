function [hfir, Ntaps] = fir_approx(Hch, fs, Efraction, N, Mct)
%% Calculate FIR approximation of channel given by anonymous function Hch
% - Inputs:
% - Hch: anonymous function of the channel frequency response
% - fs: sampling frequency of FIR filter
% - Efraction: fraction of energy that must be contained in samples of the
% channel impulse response
% - N (optional, default=1024): number of samples to use in calculations
% - Mct (optional, default=5): oversampling ratio to emulate continuous
% time

if not(exist('N', 'var'))
    N = 1024;
end

if not(exist('Mct', 'var'))
    Mct = 5;
end

f = freq_time(N, Mct*fs);
ht = fftshift(ifft(ifftshift(Hch(f))));
t0 = N/2+1;
hsamp_pos = ht(t0:Mct:end);
hsamp_neg = ht(t0:-Mct:1);
Sp = cumsum(abs(hsamp_pos).^2)/sum(abs(hsamp_pos).^2);
Sn = cumsum(abs(hsamp_neg).^2)/sum(abs(hsamp_neg).^2);
idxp = find(Sp >= Efraction, 1);
idxn = find(Sn >= Efraction, 1);
Ntaps = idxp + idxn + 1;
hfir = ht(t0 + Mct*(-idxn:idxp));