%% Validate clipping equations
clear, clc, close all

varxc = @(r) (1 + r.^2.*qfunc(r)).*(1 - qfunc(r)) - r/sqrt(2*pi).*exp(-r.^2/2).*(1 - 2*qfunc(r)) - 1/(2*pi)*exp(-r.^2);
meanxc = @(r) r.*(1 - qfunc(r)) + 1/sqrt(2*pi)*exp(-r.^2/2);

clip_ratio = 1:0.1:4;
sig = 10;
x = sig*randn(1, 2^20);

for k = 1:length(clip_ratio)
    r = clip_ratio(k);
    
    xc = x + r*sig;
    xc(xc < 0) = 0;
    
    uxc(k) = mean(xc);
    exc2(k) = mean(abs(xc).^2);
    vxc(k) = var(xc);
end

figure
plot(clip_ratio, uxc, clip_ratio, sig*meanxc(clip_ratio))

figure
plot(clip_ratio, vxc, clip_ratio, sig^2*varxc(clip_ratio))
