%% Validate expressions for the 1st and 2nd moments of a clipped normal-distributed signal

N = 2^20;       % number of points
s = 2.5476;     % std
r = 1;           % clipping ratio



for k = 1:100
    rng('shuffle'); 
    x = randn(1, N);
    xc = x;
    xc(x < -r*s) = -r*s;
    xc = xc + r*s;
    Exc(k) = mean(xc);
end

figure
stem(Exc, 'fill')

[mean(Exc), s*(r*(1-qfunc(r)) + 1/sqrt(2*pi)*exp(-r^2/2))]



   


