%% Calculate pdf using saddlepoint approximation.
% px = pdf of the output decision variable at x for all transmitted symbols
% xn = signal component projects onto orthornomal basis phi_n(t)
% varASE = ASE noise variance
% varTher = thermal noise variance
% sim = simulation struct
function px = tail_saddlepoint_approx(x, D, xn, varASE, varTher, tail)   
    [shat, ~, exitflag] = fzero(@(s) dexpoent(-abs(s), x, D, xn, varASE, varTher), 1e-3);

    if exitflag ~= 1
        warning('(x = %g) resulted in exitflag = %d\n', x, exitflag);
    end

    Ksx = expoent(-abs(shat), x, D, xn, varASE, varTher);
    d2Ksx = ddexpoent(-abs(shat), D, xn, varASE, varTher);

    px = exp(Ksx)/sqrt(2*pi*d2Ksx);
    
    if strcmp(tail, 'right')
        px = 1 - px;
    end      
    
end

% Expoent, and its first and second derivatives 
function Ksx = expoent(s, x, D, xnt, varASE, varTher)
    Ksx = sum(-log(1-D*varASE*s) + (D.*abs(xnt).^2*s)./(1-D*varASE*s)) + 0.5*varTher*s^2 - s*x - log(abs(s));
end
    
function d1Ksx = dexpoent(s, x, D, xnt, varASE, varTher) 
    d1Ksx = sum(D.*(abs(xnt).^2 + varASE - varASE^2*s)./((1-D*varASE*s).^2)) + varTher*s - x -1/s;
end

function d2Ksx = ddexpoent(s, D, xnt, varASE, varTher) 
    d2Ksx =  sum((D*varASE).^2./(1 - D*varASE*s).^2 + (2*varASE*(D.^2).*abs(xnt).^2)./((1 - D*varASE*s).^3)) + varTher +1/s^2;
end
    