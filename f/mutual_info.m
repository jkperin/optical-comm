function MI = mutual_info(X, Y, Nbins)
%% Compute mutual information between input X and output Y for a channel 
%% with conditional pdf p(y|x) when the input X is uniformly distributed 
%% over an alphabet C
% Inputs:
% - X: input to the channel. X is real-valued N x Ndim matrix, where N is 
% the number of samples and Ndim is the number of dimensions
% - Y: output to the channel. Y is real-valued N x Ndim matrix, where N is 
% the number of samples and Ndim is the number of dimensions
% - Nbins (optional, default=30): number of bins used in histogram
% calculation

% Number of bins used in histogram 
if not(exist('Nbins', 'var'))
    Nbins = 30;
end

% Compute matrix Ym. Ym is a N(m) x Ndim matrix
[C, ia, ic] = unique(X, 'rows');
[M, Ndim] = size(C);

Ym = cell(1, M);
for m = 1:M % for each symbol
    Ym{m} = Y(ic == m, :); % N(m) x Ndim matrix for the mth symbol
    N(m) = size(Ym{m}, 1);
end

% Compute edges of histogram in each dimension
dy = 1;
for n = 1:Ndim   
    ymin = min(Y(:, n));
    ymax = max(Y(:, n));
    edges{n} = linspace(ymin, 1.01*ymax, Nbins+1); % 1.01*ymax is due to a bug in function histcn
    dy = dy*abs(edges{n}(1)-edges{n}(2));
end

% Compute conditional joint pdf
pyx = cell(1, M);
lpyx = cell(1, M);
spyx = 0;
for m = 1:M % for each symbol
%     count = histcn(Ym{m}, edges{:}); % Matlab version (Mathworks contrib) 
    count = ndhistc(Ym{m}, edges{:}); % C version (Mathworks contrib) 
    pyx{m} = count/(N(m)*dy);
    lpyx{m} = log2(pyx{m} + realmin); % realmin avoids -Inf in logarithm
    spyx = spyx + pyx{m};
end
lspyx = log2(spyx + realmin); % realmin avoids -Inf in logarithm

MI = log2(M);
for m= 1:M
    MI = MI - 1/M*sum(pyx{m}(:).*(lspyx(:) - lpyx{m}(:)))*dy;
end



% 
% 
% 
% 
% if not(all(isreal(Y)))
%     Yr = real(Y);
%     Yi = imag(Y);
% 
%     
%     yimin = min(min(Yi));
%     yimax = max(max(Yi));
% 
%     
% 
% 
%     pyx = cell(1, M);
%     lpyx = cell(1, M);
%     spyx = 0;
%     for k = 1:M
% %         pyx{k} = histcounts2(Yr(:, k), Yi(:, k), redges, iedges, 'Normalization', 'pdf');
%         count = histcn([Yr(:, k), Yi(:, k)], redges, iedges);
%         pyx{k} = count/(N*dyr*dyi);
%         lpyx{k} = log2(pyx{k} + realmin); % realmin avoids -Inf in logarithm
%         spyx = spyx + pyx{k};
%     end
%     lspyx = log2(spyx + realmin); % realmin avoids -Inf in logarithm
% 
%     MI = log2(M);
%     for k= 1:M
%         MI = MI - 1/M*sum(sum(pyx{k}.*(lspyx - lpyx{k})))*dyr*dyi;
%     end
% else
%     ymin = min(min(Y));
%     ymax = max(max(Y));
% 
%     edges = linspace(ymin, ymax, Nbins+1);
% 
%     pyx = cell(1, M);
%     lpyx = cell(1, M);
%     spyx = 0;
%     for k = 1:M
%         pyx{k} = histcounts(Y(:, k), edges, 'Normalization', 'pdf');
%         lpyx{k} = log2(pyx{k} + realmin); % realmin avoids -Inf in logarithm
%         spyx = spyx + pyx{k};
%     end
%     lspyx = log2(spyx + realmin); % realmin avoids -Inf in logarithm
% 
%     dy = abs(edges(1)-edges(2));
% 
%     MI = log2(M);
%     for k= 1:M
%         MI = MI - 1/M*sum(pyx{k}.*(lspyx - lpyx{k}))*dy;
%     end
% end
%     