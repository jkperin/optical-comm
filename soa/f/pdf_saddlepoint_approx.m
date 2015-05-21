%% Calculate pdf using saddlepoint approximation.
% px = pdf of the output decision variable at x for all transmitted symbols
% xn = signal component projects onto orthornomal basis phi_n(t)
% varASE = ASE noise variance
% varTher = thermal noise variance
% sim = simulation struct
function px = pdf_saddlepoint_approx(x, D, xn, varASE, varTher, sim)
    if sim.verbose
        figure, hold on
    end
    
    px = zeros(length(x), sim.Nsymb);
    for k = 1:sim.Nsymb
        for kk = 1:length(x)
            [shat, ~, exitflag] = fzero(@(s) dexpoent(s, x(kk), D, xn(k, :).', varASE, varTher), 1e-3);

            if exitflag ~= 1
                warning('(%d, %g) resulted in exitflag = %d\n', k, x(kk), exitflag);
                dataTX(k)
            end

            Ksx = expoent(shat, x(kk), D, xn(k, :).', varASE, varTher);
            d2Ksx = ddexpoent(shat, D, xn(k, :).', varASE, varTher);
            
%             if d2Ksx < 0
%                 warning('(%d, %g) resulted in kp2 < 0\n', k, x(kk));
%                 d2Ksx
%             end

            px(kk, k) = real(exp(Ksx)/sqrt(2*pi*d2Ksx)); 

        end
    
        if sim.verbose
            plot(x, px(:, k))
        end
        
        px(:, k) =  px(:, k)/trapz(x,  px(:, k));
    end
end

% Expoent, and its first and second derivatives 
function Ksx = expoent(s, x, D, xnt, varASE, varTher) 
% if symb ~= 0
    Ksx = sum(-log(1-D*varASE*s) + (D.*abs(xnt).^2*s)./(1-D*varASE*s)) + 0.5*varTher*s^2 - s*x;
% else
%     Ksx = sum(log(1./(1-D*varASE*s))) + 0.5*varTher*s^2 - s*x;
% end
end
    
function d1Ksx = dexpoent(s, x, D, xnt, varASE, varTher) 
% if symb ~= 0
    d1Ksx = sum(D.*(abs(xnt).^2 + varASE - varASE^2*s)./((1-D*varASE*s).^2)) + varTher*s - x;
% else
%     d1Ksx = sum((D.*varASE)./(1-D*varASE*s)) + varTher*s - x;
% end
1;
end

function d2Ksx = ddexpoent(s, D, xnt, varASE, varTher) 
% if symb ~= 0
    d2Ksx = sum((D*varASE).^2./(1 - D*varASE*s).^2 + (2*varASE*(D.^2).*abs(xnt).^2)./((1 - D*varASE*s).^3)) + varTher;
% else
%     d2Ksx = sum((D*varASE).^2./(1 - D*varASE*s).^2) + varTher;
% end
1;
end
    