function px = pdf_saddlepoint_approx(x, xn, varASE, varTher, sim)

    figure, hold on
    px = zeros(sim.Npoints, sim.Me);
    for k = 1:sim.Nsymb
        for kk = 1:length(x)
            [shat, ~, exitflag] = fzero(@(s) dexpoent(s, x(kk), D, xn(k, :).', varASE, varTher), 0);

            if exitflag ~= 1
                warning('(%d, %g) resulted in exitflag = %d\n', k, x(kk), exitflag);
                dataTX(k)
            end

            Ksx = expoent(shat, x(kk), D, xnd(k, :).', varASE, varTher);
            d2Ksx = ddexpoent(shat, D, xnd(k, :).', varASE, varTher);
            
            if d2Ksx > 0
                warning('(%d, %g) resulted in kp2 < 0\n', k, x(kk));
                d2Ksx
            end

            px(kk, k) = exp(Ksx)/sqrt(-2*pi*d2Ksx); 

        end
    
        if dataTX(k) ~= 0
            plot(x, px(:, k))
        end
    end
end

% Expoent, and its first and second derivatives 
function Ksx = expoent(s, x, D, xnt, varASE, varTher) 
    Ksx = sum(log(1./(1-D*varASE*s)) + (D.*abs(xnt).^2*s)./(1-D*varASE*s)) + 0.5*varTher*s^2 - s*x;
    1;
end
    
function d1Ksx = dexpoent(s, x, D, xnt, varASE, varTher) 
    d1Ksx = sum(D.*(abs(xnt).^2 + varASE - varASE^2*s)./((1-D*varASE*s).^2)) + varTher*s - x;
    1;
end

function d2Ksx = ddexpoent(s, D, xnt, varASE, varTher) 
    d2Ksx = sum((D*varASE).^2./(1 - D*varASE*s).^2 - (2*varASE*(D.^2).*abs(xnt).^2)./((1 - D*varASE*s).^3)) + ...
        varTher;
    1;
end
    