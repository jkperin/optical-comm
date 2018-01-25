function D = GN_model_coeff_qsub(spanLengthkm, l_idx)
%% Compute GN model coefficients on cluster

if any(not(isnumeric([spanLengthkm l_idx])))
    spanLength = 1e3*str2double(spanLengthkm);
    l_idx = round(str2double(l_idx));
end

Fiber = fiber(spanLength, @(l) 0.18);

Df = 33.3e9;
dlamb = df2dlamb(Df);
lamb = 1525e-9:dlamb:1570e-9;

% Df = 50e9;
% dlamb = df2dlamb(Df);
% lamb = 1522e-9:dlamb:1582e-9;

% lcenter = 1550e-9;
% fcenter = lambda2freq(lcenter);
% f = (-200:50:200)*1e9 + fcenter;
% lamb = freq2lambda(f);

%% Test GN model
tic
D = GN_model_coeff(lamb, Df, Fiber, l_idx);
toc

filename = sprintf('GN_model_coeff_spanLengthkm=%s_Df=%dGHz_l=%d.mat', spanLengthkm, round(Df/1e9), l_idx);
save(filename)

% figure, imagesc(10*log10(D))
