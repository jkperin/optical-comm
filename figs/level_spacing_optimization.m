%% Level spacing optimization
clear, clc, close all
M = 16; % number of levels
rexdB = -Inf; % extinction ratio
var1 = 1; % variance of signal-independ noise
vars_ratio = 0.5; % 
BERtarget = 1e-4;

% [a, b] = level_optimization1(M, BERtarget, rexdB, var1, vars_ratio);
% 
% figure(1), hold on
% x = linspace(a(1)-5, a(end)+10);
% var2 = @(Plevel) (var1*vars_ratio)*Plevel;
% sigk = @(Plevel) sqrt(var1 + var2(Plevel));
% for level = 1:M-1
%     plot(x, pdf('normal', x, a(level), sigk(a(level))))
%     plot(x, pdf('normal', x, a(level+1), sigk(a(level+1))))
%     plot(b(level)*[1 1], [0 0.1], 'k')
% end

r = logspace(0, 2, 20);
A1 = zeros(size(r));
A2 = zeros(size(r));
for k = 1:length(r)
    vars_ratio = r(k); % sigma_{dependent}/(sigma_{independent}P)
    [a1, ~] = level_optimization1(M, BERtarget, rexdB, var1, vars_ratio);
    [a2, ~] = level_optimization2(M, BERtarget, rexdB, var1, vars_ratio);
    
    A1(k) = mean(a1);
    A2(k) = mean(a2);    
end

figure, hold on, box on
plot(r, A1)
plot(r, A2)
legend('Algorithm 1', 'Algorithm 2', 'Location', 'NorthWest')
xlabel('\sigma^2_{dependent}/(\sigma^2_{independent}\cdot P)')

figure, box on
plot(r, 10*log10(A1./A2))
xlabel('\sigma^2_{dependent}/(\sigma^2_{independent}\cdot P)')
ylabel('Average power improvement (dB)')

    

    

    
% figure, hold on
% plot(log(tol))
% plot([1 k], log(maxtol*[1 1]), 'r')
% xlabel('Iteration')
% ylabel('log(Tolerance)')
% legend('Tolerance', 'Required for Convergence')
% title('Level optimization convergece')