function [a, b] = level_spacing_optm(mpam, tx, soa, rx, sim)

a = zeros(mpam.M, 1);
b = zeros(mpam.M-1, 1);

a(1) = tx.Pmin*soa.Gain;

Pe = log2(mpam.M)*sim.BERtarget*mpam.M; % error probability given a level

% KLSE Fourier Series Expansion 
[U, D, Fmax] = klse_fourier(rx, sim, sim.Mct); 

df = 1/sim.Mct;
ft = (-0.5:df:0.5-df).';
f = ft(abs(ft) <= Fmax);
fm = f(abs(f) <= Fmax);

% Noises variance
varASE = (soa.N0*sim.fs/2)/sim.Mct; % variance of circularly complex Gaussain noise 
varTherm = rx.N0*rx.elefilt.noisebw(sim.fs)/2;

for level = 1:mpam.M-1
    % Find decision threshold
    ck = calc_mgf_param(a(level), U, fm);
    
    [dPthresh, ~, exitflag] = fzero(@(dPthresh)...
        tail_saddlepoint_approx(a(level) + abs(dPthresh), D, ck, varASE, varTherm, 'right') - Pe,...
        sqrt(varASE)*qfuncinv(Pe));
    
    if exitflag ~= 1
        warning('level_spacing_optm: threshold optimization did not converge');
    end
    
    b(level) = a(level) + abs(dPthresh);
    
    % Find next level  
    [dPlevel, ~, exitflag] = fzero(@(dPlevel)...
        tail_saddlepoint_approx(b(1), D, calc_mgf_param(b(level) + abs(dPlevel), U, fm), varASE, varTherm, 'left') - Pe,...
        2*b(level)-a(level));    
    
    if exitflag ~= 1
        warning('level_spacing_optm: level optimization did not converge');     
    end
    
    a(level+1) = b(level) + abs(dPlevel);
end

function ck = calc_mgf_param(Plevel, U, fm)
% Calculate output
xk = zeros(size(fm));
xk(fm == 0) = sqrt(Plevel);

% Used to calculate non-centrality parameter of Chi-Squared distributions
% ck(i, j) = ith chi-square distribution at jth time sample. |ck(i, j)|^2
% is the non-centrality parameter of that chi-square distribution
ck = U'*xk; 


