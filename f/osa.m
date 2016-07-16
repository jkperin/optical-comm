function [YdBm, ly, fy] = osa(x, lambda, f, resolutionnm, verbose)
%% Optical spectrum analyser
% Inputs:
% - x: input signal
% - lambda: wavelength of the input signal
% - f: frequency vector
% - resolution (optional, default = 0.1 nm): resolution in nm
% - verbose (optional, default = false): whether to plot spectrum at the
% end

c = 299792458; % speed of light
if exist('resolutionnm', 'var') && resolutionnm
    resolution = resolutionnm*c/((lambda+resolutionnm/2)*(lambda-resolutionnm/2));
else
    resolutionnm = 0.1e-9;
    resolution = 12.5e9; % default 0.1nm at 1550nm
end
    
if size(x, 1) == 1 % 1 pol
    P = abs(fftshift(fft(x))).^2;
else % 2-pol
    P = abs(fftshift(fft(x(1, :)))).^2 + abs(fftshift(fft(x(2, :)))).^2;
end

fs = -2*f(1);
P = P/(fs*length(P)/resolution);

df = abs(f(2) - f(1));
windowLength = round(resolution/df);
if mod(windowLength, 2) == 0 % if even
    windowLength = windowLength + 1;
end

Pfilt = filter(ones(1, windowLength)/windowLength, 1, P);
Pfilt = circshift(Pfilt, [0 -(windowLength-1)/2]); % remove delay due to filter

fc = c/lambda;
fy = fc + [fliplr(0:-resolution:min(f)) resolution:resolution:max(f)];
Y = interp1(f, Pfilt, fy-fc);
ly = c./fy;
YdBm = 10*log10(Y/1e-3);

if exist('verbose', 'var') && verbose
    figure(132), box on
    plot(ly*1e9, YdBm)
    xlabel('Wavelength (nm)')
    ylabel('Power (dBm)')
    title(sprintf('Optical spectrum with resolution %.2f nm', 1e9*resolutionnm))
end    


