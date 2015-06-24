%% Level spacing (a) and decision threshold (b) optmization
% Assumes infinite extinction ratio at first, then corrects power and
% optmize levels again
% The levels and thresholds calculated are after APD amplification
function [a, b] = level_spacing_optm_gauss_approx(mpam, tx, apd, rx, sim)

% Error probability under a single tail for a given symbol
Pe = log2(mpam.M)*sim.BERtarget*(mpam.M/(2*(mpam.M-1)));

% Noises variance
varTherm = rx.N0*rx.elefilt.noisebw(sim.fs)/2; % variance of thermal noise
% Shot noise variance = Agrawal 4.4.17 (4th edition)
calc_noise_std = @(Plevel) sqrt(varTherm + 2*apd.q*apd.Gain^2*apd.Fa*(apd.R*Plevel/apd.Gain + apd.Id)*rx.elefilt.noisebw(sim.fs)/2);

% Initialize levels and thresholds
a = zeros(mpam.M, 1);
b = zeros(mpam.M-1, 1);

rex = 10^(-abs(tx.rexdB)/10);

maxtol = 1e-6; % maximum tolerance for convergence
maxit = 20; % maximum number of iteratios

tol = Inf;
k = 1;
while tol(end) > maxtol && k < maxit

    apast = a;
    a(1) = a(end)*rex;

    for level = 1:mpam.M-1
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

if sim.verbose
    figure, hold on
    plot(log(tol))
    plot([1 k], log(maxtol*[1 1]), 'r')
    xlabel('Iteration')
    ylabel('log(Tolerance)')
    legend('Tolerance', 'Required for Convergence')
    title('Level optimization convergece')
end