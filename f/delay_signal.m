function y = delay_signal(x, delay)
%% Delay signal by delay samples, which need not be integer

transpose = false;
if size(x, 1) > size(x, 2)
    x = x.';
    transpose = true;
end

if isInteger(delay) % just circshift signal
    delay = round(delay);
    y = circshift(x, [0 delay]);
else
    df = 1/length(x);
    f = ifftshift(-0.5:df:0.5-df);    
    y = ifft(fft(x).*exp(-1j*2*pi*f*delay));
end

if all(isreal(x))
    y = real(y);
end

if transpose
    y = y.';
end
    