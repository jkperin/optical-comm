%% Validate saddlepoint approximation with Gaussian
clear, clc, close all
addpath f

mu = 2;

mu1 = 2;
mu2 = 5;
sig2 = 5;

x = linspace(0, 100, 100);

for k= 1:length(x)
%       Gaussian
%     [shat(k),fval, exitflag] = fzero(@(s) mu + sig2*s - x(k), 0);
%     shat2(k) = (x(k)-mu)/sig2;
% 
%     Ks = mu*shat2(k) + 1/2*sig2*shat2(k)^2-shat2(k)*x(k);
%     Ks2 = sig2;
%  
%     Non-central Chi-square
%     [shat(k),fval, exitflag] = fzero(@(s) sig2/(1-sig2*s) + (mu1^2+mu2^2)/(1-sig2*s)^2 - x(k), 0);
    varTher = 0;
    D = 1;
    xnt = mu1+1j*mu2;
    varASE = sig2;
    [shat(k),fval, exitflag] = fzero(@(s) sum(D.*(abs(xnt).^2 + varASE - varASE^2*s)./((1-D*varASE*s).^2)) + varTher*s - x(k), 0);
    
    
%     Ks(k) = -log(1-sig2*shat(k)) + (mu1^2+mu2^2)*shat(k)/(1-sig2*shat(k)) -x(k)*shat(k);
    Ks(k) = sum(-log(1-D*varASE*shat(k)) + (D.*abs(xnt).^2*shat(k))./(1-D*varASE*shat(k))) + 0.5*varTher*shat(k)^2 - shat(k)*x(k);
    
    Ks2(k) = sum((D*varASE).^2./(1 - D*varASE*shat(k)).^2 + (2*varASE*(D.^2).*abs(xnt).^2)./((1 - D*varASE*shat(k)).^3)) + varTher;
%     Ks2(k) = sig2^2/(1-sig2*shat(k))^2 + 2*sig2*(mu1^2+mu2^2)/(1-sig2*shat(k))^3;
    
    dKsx = @(s) sum(D.*(abs(xnt).^2 + varASE - varASE^2*s)./((1-D*varASE*s).^2)) + varTher*s - x(k) -1/s;
    Ksx =@(s) sum(-log(1-D*varASE*s) + (D.*abs(xnt).^2*s)./(1-D*varASE*s)) + 0.5*varTher*s^2 - s*x(k) - log(abs(s));
    ddKsx = @(s) sum((D*varASE).^2./(1 - D*varASE*s).^2 + (2*varASE*(D.^2).*abs(xnt).^2)./((1 - D*varASE*s).^3)) + varTher +1/s^2;

%     Tail probabilities
    [shatt(k),fval, exitflag] = fzero(@(s) dKsx(-abs(s)), 0.001);
    
%     Ks(k) = -log(1-sig2*shat(k)) + (mu1^2+mu2^2)*shat(k)/(1-sig2*shat(k)) -x(k)*shat(k);
    Kstail(k) = Ksx(-abs(shatt(k)));
    
    Ks2tail(k) = ddKsx(-abs(shatt(k)));
%     Ks2(k) = sig2^2/(1-sig2*shat(k))^2 + 2*sig2*(mu1^2+mu2^2)/(1-sig2*shat(k))^3;
    
    % Saddlepoint approximation
    px(k) = exp(Ks(k))/sqrt(2*pi*Ks2(k));
    
    pxtail(k) = exp(Kstail(k))/sqrt(2*pi*Ks2tail(k));
      
end

x1 = sqrt(sig2/2)*randn(2^14, 1) + mu1;
x2 = sqrt(sig2/2)*randn(2^14, 1) + mu2;

z = abs(x1).^2 + abs(x2).^2;

% ptrue = pdf('Normal', x, mu, sqrt(sig2));
ptrue = 1/(sig2/2)*pdf('Noncentral Chi-square', x/(sig2/2), 2, (mu1^2+mu2^2)/(sig2/2));
cdftrue = cdf('Noncentral Chi-square', x/(sig2/2), 2, (mu1^2+mu2^2)/(sig2/2));

sim.Nsymb = 1;
sim.verbose = false;
ptest = pdf_saddlepoint_approx(x, 1, xnt, sig2, 0, sim);
% ptest = ptest/trapz(x, ptest);
for k = 1:length(x)
    ptail_test(k) = tail_saddlepoint_approx(x(k), D, xnt, varASE, varTher, 'left');  
end

figure
% plot(x, shat, x, shat2, '--')
 plot(x, shat)
figure, hold on
plot(x, px, x, ptrue)
plot(x, ptest, '--k')
% [n,xb] = hist(z, 50);
% n = n/trapz(xb, n);
% bar(xb, n)
legend('saddle-point approx', 'theory', 'program', 'count')

figure, hold on
plot(x, pxtail, '.-k', x, ptail_test)
plot(x, cdftrue)
% 
% figure
% plot(x, shatt)
