% Given signal x, compute the spectrum Sxx using a periodogram with block size Nb

function [X, f] = spectrum_periodogram(x, Nfft, fs)

N_symb = length(x)/Nfft; % number of symbols to average

x = reshape(x(1:N_symb*Nfft), N_symb, Nfft);
X = fftshift(fft(x, Nfft, 1))/Nfft;
X = mean(abs(X).^2, 2);

df = fs/Nfft;
f = -fs/2:df:fs/2-df;

plot(f/1e9, 10*log10(X))
xlabel('Frequency (GHz)')
ylabel('Spectrum from periodogram')
a = axis;
axis([-fs/2e9 fs/2e9 a(3) a(4)]);




