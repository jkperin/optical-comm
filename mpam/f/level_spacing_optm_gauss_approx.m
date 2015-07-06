%% Level spacing (a) and decision threshold (b) optmization
% Assumes infinite extinction ratio at first, then corrects power and
% optmize levels again
% The levels and thresholds calculated are after APD amplification
function [a, b] = level_spacing_optm_gauss_approx(M, BERtarget, rexdB, calc_noise_std, verbose)

% Error probability under a single tail for a given symbol
Pe = log2(M)*BERtarget*(M/(2*(M-1)));

% Initialize levels and thresholds
a = zeros(M, 1);
b = zeros(M-1, 1);

rex = 10^(-abs(rexdB)/10);

maxtol = 1e-6; % maximum tolerance for convergence
maxit = 20; % maximum number of iteratios

tol = Inf;
k = 1;
while tol(end) > maxtol && k < maxit

    apast = a;
    a(1) = a(end)*rex;

    for level = 1:M-1
        % Find threshold
        sig = calc_noise_std(a(level));
        
        [dPthresh, ~, exitflag] = fzero(@(dPthresh) qfunc(abs(dPthresh)/sig) - Pe, 0);

        if exitflag ~= 1
            warning('level_spacing_optm: threshold optimization did not converge');
        end

        b(level) = a(level) + abs(dPthresh);

        % Find next level  
        [dPlevel, ~, exitflag] = fzero(@(dPlevel) qfunc(abs(dPlevel)/calc_noise_std(b(level) + abs(dPlevel))) - Pe, 0);    

        if exitflag ~= 1
            warning('level_spacing_optm: level optimization did not converge');     
        end

        a(level+1) = b(level) + abs(dPlevel);
    end
    
    tol(k) = sqrt(sum(abs(a-apast).^2));
    k = k + 1;       
end

if nargin == 5 && verbose
    figure, hold on
    plot(log(tol))
    plot([1 k], log(maxtol*[1 1]), 'r')
    xlabel('Iteration')
    ylabel('log(Tolerance)')
    legend('Tolerance', 'Required for Convergence')
    title('Level optimization convergece')
end