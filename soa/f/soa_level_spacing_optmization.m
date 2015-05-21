function [a, b] = soa_level_spacing_optmization(mpam, tx, soa, rx, OptFilt, EleFilt, sim)

a(1) = tx.Pmin;

Pe = log2(mpam.M)*sim.BERtarget*mpam.M;

for k = 1:mpam.M-1
    [px, x] = px_soa_klse_freq(a(k), soa, rx, OptFilt, EleFilt, sim);
    1;
    
    [xth, fval, exitflag] = fzero(@(xth) trapz(x(x > soa.Gain*(a(k) + abs(xth))), px(x > soa.Gain*(a(k) + abs(xth)))) - Pe, 0);
    
    if exitflag ~= 1
        warning('did not converge');
    end
    
    b(k) = a(k) + abs(xth);
    
    [xi, fval, exitflag] = fzero(@(xi) calc_ptail_fix_thresh(b(k) + abs(xi), soa.Gain*b(k), soa, rx, OptFilt, EleFilt, sim) - Pe, 0);

    a(k+1) = b(k) + abs(xi);
end