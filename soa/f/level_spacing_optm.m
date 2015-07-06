%% Level spacing (a) and decision threshold (b) optmization
% Assumes infinite extinction ratio at first, then corrects power and
% optmize levels again
% The calculated levels and thresholds are at the receiver
function [a, b] = level_spacing_optm(mpam, tx, soa, rx, sim)

N = sim.Mct+1;

% Error probability under a single tail for a given symbol
Pe = log2(mpam.M)*sim.BERtarget*(mpam.M/(2*(mpam.M-1)));

% KLSE Fourier Series Expansion 
[U, D, Fmax] = klse_fourier(rx, sim, N); 

df = 1/N;
ft = (-0.5:df:0.5-df).';
f = ft(abs(ft) <= Fmax);
fm = f(abs(f) <= Fmax);

% Noises variance
varASE = (soa.N0*sim.fs)/N; % variance of circularly complex Gaussain noise 
% Note: even though soa.N0 is single-sided PSD we don't multiply by
% sim.fs/2 because this is a band-pass process
varTherm = rx.N0*rx.elefilt.noisebw(sim.fs)/2;

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
        % Find decision threshold
        ck = calc_mgf_param(a(level), U, fm);

        try
        [dPthresh, ~, exitflag] = fzero(@(dPthresh)...  
            tail_saddlepoint_approx(a(level) + abs(dPthresh), D, ck, varASE, varTherm, 'right') - Pe, 0);
        catch e
            disp(e.message)
            1;
        end

        if exitflag ~= 1
            warning('level_spacing_optm: threshold optimization did not converge');
        end

        b(level) = a(level) + abs(dPthresh);

        % Find next level  
        [dPlevel, ~, exitflag] = fzero(@(dPlevel)...
            tail_saddlepoint_approx(b(level), D, calc_mgf_param(b(level) + abs(dPlevel), U, fm), varASE, varTherm, 'left') - Pe, 0);    

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


function ck = calc_mgf_param(Plevel, U, fm)
% Calculate output
xk = zeros(size(fm));
xk(fm == 0) = sqrt(Plevel);

% Used to calculate non-centrality parameter of Chi-Squared distributions
% ck(i, j) = ith chi-square distribution at jth time sample. |ck(i, j)|^2
% is the non-centrality parameter of that chi-square distribution
ck = U'*xk; 


