function [PdBm, PW] = power_meter(x)
if size(x, 1) > size(x, 2)
    x = x.';
end

if size(x, 1) == 1 % single pol
    PW = mean(abs(x).^2);
else % dual pol
    PW = mean(abs(x(1, :)).^2) + mean(abs(x(2, :)).^2);
end

PdBm = Watt2dBm(PW);