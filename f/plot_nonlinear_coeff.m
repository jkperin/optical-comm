clear, clc, close all

%% SMF
load('GN_coeff/legacy/GN_model_coeff_spanLengthkm=50km_Df=50GHz.mat')
for n = 1:3
    nonlinear_coeff{n} = (1.4/0.8)^2*nonlinear_coeff{n};
end
clims = [-100 10*log10(max(nonlinear_coeff{2}(:)))];

figure, imagesc(10*log10(nonlinear_coeff{1}), clims)
axis square
colormap hot
set(gca, 'xtick', [1 150 299])
set(gca, 'xticklabels', {'-N+1', 0, 'N-1'})
set(gca, 'ytick', [1 150 299])
set(gca, 'yticklabels', {'N-1', 0, '-N+1'})

figure, imagesc(10*log10(nonlinear_coeff{2}), clims)
axis square
colormap hot
set(gca, 'xtick', [1 150 299])
set(gca, 'xticklabels', {'-N+1', 0, 'N-1'})
set(gca, 'ytick', [1 150 299])
set(gca, 'yticklabels', {'N-1', 0, '-N+1'})

figure, imagesc(10*log10(nonlinear_coeff{3}), clims)
axis square
colormap hot
set(gca, 'xtick', [1 150 299])
set(gca, 'xticklabels', {'-N+1', 0, 'N-1'})
set(gca, 'ytick', [1 150 299])
set(gca, 'yticklabels', {'N-1', 0, '-N+1'})
colorbar

%% Large effective area fiber
load('GN_model_coeff_spanLengthkm=50km_Df=50GHz.mat')

% clims = [-100 10*log10(max(nonlinear_coeff{2}(:)))];

figure, imagesc(10*log10(nonlinear_coeff{1}), clims)
axis square
colormap hot
set(gca, 'xtick', [1 150 299])
set(gca, 'xticklabels', {'-N+1', 0, 'N-1'})
set(gca, 'ytick', [1 150 299])
set(gca, 'yticklabels', {'N-1', 0, '-N+1'})

figure, imagesc(10*log10(nonlinear_coeff{2}), clims)
axis square
colormap hot
set(gca, 'xtick', [1 150 299])
set(gca, 'xticklabels', {'-N+1', 0, 'N-1'})
set(gca, 'ytick', [1 150 299])
set(gca, 'yticklabels', {'N-1', 0, '-N+1'})

figure, imagesc(10*log10(nonlinear_coeff{3}), clims)
axis square
colormap hot
set(gca, 'xtick', [1 150 299])
set(gca, 'xticklabels', {'-N+1', 0, 'N-1'})
set(gca, 'ytick', [1 150 299])
set(gca, 'yticklabels', {'N-1', 0, '-N+1'})
colorbar