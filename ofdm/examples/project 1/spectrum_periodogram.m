% Given signal x, compute the spectrum Sxx using a periodogram with block size Nb

function Sxx = spectrum_periodogram(x,Nb)

N = length(x);
M = floor(N/Nb);

Sxx = zeros(1,Nb);

for m = 0:M-1
    Sxx = Sxx + (1/M)*abs(fft(x(m*Nb+[1:Nb]))).^2;
end;