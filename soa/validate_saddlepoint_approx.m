%% Validate saddlepoint approximation with Gaussian and noncentral Chi-Square
clear, clc, close all
addpath f

mu = 2;

mu1 = 2;
mu2 = 5;
sig2 = 5;

x = linspace(0, 100, 100);
shat_pdf = zeros(size(x));
shat_tail =  zeros(size(x));
for k = 1:length(x)
%       Gaussian
%     [shat(k),fval, exitflag] = fzero(@(s) mu + sig2*s - x(k), 0);
%     shat2(k) = (x(k)-mu)/sig2;
% 
%     Ks = mu*shat2(k) + 1/2*sig2*shat2(k)^2-shat2(k)*x(k);
%     Ks2 = sig2;
%  
%     Non-central Chi-square
%     [shat(k),fval, exitflag] = fzero(@(s) sig2/(1-sig2*s) + (mu1^2+mu2^2)/(1-sig2*s)^2 - x(k), 0);

    % Put in same notation as SOA code
    varTher = 0;        % thermal noise variance
    D = 1;              % Number of noncentral chi square (number of eigenvalues in KL expansion)
    xnt = mu1+1j*mu2;   % noncentrality parameter written in complex form
    varASE = sig2;      % variance of two Gaussians that given orgin to chi-square
      
    %% PDF (no renormalization is done after saddlepoint approximation)
    % K(s, x) and its second derivatives for saddlepoint approximation of
    % pdf. In this case K(s, x) = log(M(s)) - sx.
    dKsx = @(s, x) sum(D.*(abs(xnt).^2 + varASE - varASE^2*s)./((1-D*varASE*s).^2)) + varTher*s - x;
    Ksx = @(s, x) sum(-log(1-D*varASE*s) + (D.*abs(xnt).^2*s)./(1-D*varASE*s)) + 0.5*varTher*s^2 - s*x;
    ddKsx = @(s) sum((D*varASE).^2./(1 - D*varASE*s).^2 + (2*varASE*(D.^2).*abs(xnt).^2)./((1 - D*varASE*s).^3)) + varTher;
    
    % Calculate saddlepoint shat for pmf approximation
    [shat_pdf(k), fval, exitflag] = fzero(@(s) dKsx(s, x(k)), 0);
       
    px(k) = exp(Ksx(shat_pdf(k), x(k)))/sqrt(2*pi*ddKsx(shat_pdf(k)));
    
    %% Tail
    % K(s, x) and its first two derivatives for saddlepoint approximation
    % of tail probabilities. In this case K(s, x) = log(M(s)/s) - sx
    dKsx = @(s, x) sum(D.*(abs(xnt).^2 + varASE - varASE^2*s)./((1-D*varASE*s).^2)) + varTher*s - x -1/s;
    Ksx =@(s, x) sum(-log(1-D*varASE*s) + (D.*abs(xnt).^2*s)./(1-D*varASE*s)) + 0.5*varTher*s^2 - s*x - log(abs(s));
    ddKsx = @(s) sum((D*varASE).^2./(1 - D*varASE*s).^2 + (2*varASE*(D.^2).*abs(xnt).^2)./((1 - D*varASE*s).^3)) + varTher +1/s^2;

    % Tail probabilities
    % enforce negative saddlepoint in order to calculate left tail
    % (x-> -\infty). As a result, tail will only be accurate if x << E(x)
    [shat_tail(k),fval, exitflag] = fzero(@(s) dKsx(-abs(s), x(k)), 0.001);
    shat_tail(k) = -abs(shat_tail(k));
        
    pxtail(k) = exp(Ksx(shat_tail(k), x(k)))/sqrt(2*pi*ddKsx(shat_tail(k)));

end

% Generate z distributed according to noncentral chi2
x1 = sqrt(sig2/2)*randn(2^14, 1) + mu1;
x2 = sqrt(sig2/2)*randn(2^14, 1) + mu2;

z = abs(x1).^2 + abs(x2).^2;

% True pdf and cdf
% ptrue = pdf('Normal', x, mu, sqrt(sig2));
ptrue = 1/(sig2/2)*pdf('Noncentral Chi-square', x/(sig2/2), 2, (mu1^2+mu2^2)/(sig2/2));
cdftrue = cdf('Noncentral Chi-square', x/(sig2/2), 2, (mu1^2+mu2^2)/(sig2/2));

mean_true = trapz(x, x.*ptrue);

% Generate pdf using saddlepoint approximation (with renormalization)
ptest = pdf_saddlepoint_approx(x, 1, xnt, sig2, 0);

% Generate tail probabilities using saddlepoint approximation
ptail_test = tail_saddlepoint_approx(x(x < mean_true), D, xnt, varASE, varTher, 'left'); 
ptail_test = [ptail_test 1-tail_saddlepoint_approx(x(x >= mean_true), D, xnt, varASE, varTher, 'right')]; 

ptail_pos = tail_saddlepoint_approx(x, D, xnt, varASE, varTher, 'right'); 
ptail_pos = 1 - ptail_pos;

figure
plot(x, shat_pdf, x, shat_tail)
legend('pdf', 'tail')
ylabel('saddlepoint')
xlabel('x')

figure, hold on
plot(x, px, x, ptest, '--')
plot(x, ptrue, 'k')
% [n,xb] = hist(z, 50);
% n = n/trapz(xb, n);
% bar(xb, n)
legend('saddlepoint approx', 'saddlepoint approx (w/ renormalization)', 'theory')

figure, hold on
plot(x, pxtail, 'b', x, ptail_pos, 'r', x, ptail_test, '.-m')
plot(x, cdftrue, 'k')
legend('tail with negative saddleppoint', 'tail with positive saddleppoint', 'tail alternating saddleppoint', 'true cdf', 'Location', 'SouthEast')

