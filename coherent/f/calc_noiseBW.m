function noiseBW = calc_noiseBW(Wx, Wy)
%% Calculates one-sided noise bandwidth from discrete-time filters with coefficients Wx and Wy
% NoiseBW is the average noise bandwidth of all the filters in Wx and Wy

Nfilters = size(Wx, 2);
noiseBWx = zeros(1, Nfilters);
noiseBWy = zeros(1, Nfilters);
w = 2*pi*linspace(-0.5, 0.5, 1e3);
for k = 1:size(Wx, 2)
    H = freqz(Wx(:, k), 1, w);
    noiseBWx(k) = 0.5*trapz(w/(2*pi), abs(H).^2)/abs(interp1(w, H, 0)).^2;
    
    H = freqz(Wy(:, k), 1, w);
    noiseBWy(k) = 0.5*trapz(w/(2*pi), abs(H).^2)/abs(interp1(w, H, 0)).^2;
end

noiseBW = mean([noiseBWx noiseBWy]);