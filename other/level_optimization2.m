function [a,b] = level_optimization2(M, BERtarget, rexdB, var1, vars_ratio)
    %% Algorithm 2
    Pe = log2(M)*BERtarget*M/(M-1);
    tol = Inf;
    k = 1;
    var2 = @(Plevel) (var1*vars_ratio)*Plevel;
    sigk = @(Plevel) sqrt(var1 + var2(Plevel));
    
    rex = 10^(-abs(rexdB)/10);
    
    a = zeros(M, 1);
    b = zeros(M-1, 1);

    maxtol = 1e-6; % maximum tolerance for convergence
    maxit = 20; % maximum number of iteratios
    while tol(end) > maxtol && k < maxit

        apast = a;
        a(1) = a(end)*rex;

        for level = 1:M-1
            % Find next level  

            [dP, ~, exitflag] = fzero(@(dP) qfunc((intersection_point(a(level), a(level) + abs(dP)) - a(level))/sigk(a(level)))...
                + qfunc((a(level) + abs(dP) - intersection_point(a(level), a(level) + abs(dP)))/sigk(a(level) + abs(dP))) - Pe, sigk(a(level)));    

            if exitflag ~= 1
                warning('level_spacing_optm: level optimization did not converge');     
            end

            a(level+1) = a(level) + abs(dP);
            b(level) = intersection_point(a(level), a(level+1));

%             figure(1), hold on
%             x = linspace(a(1)-5, a(end)+100);
%             plot(x, pdf('normal', x, a(level), sigk(a(level))))
%             plot(x, pdf('normal', x, a(level+1), sigk(a(level+1))))
%             plot(b(level)*[1 1], [0 0.1], 'k')
%             1;
        end

        tol(k) = sqrt(sum(abs(a-apast).^2));
        k = k + 1;  
    end

    Nsymb = 2^20;
    dataTX = randi([0 M-1], 1, Nsymb);

    x = a((gray2bin(dataTX, 'pam', M) + 1));

    y = x + sqrt(var1 + var2(x)).*randn(Nsymb, 1);

    dataRX = sum(bsxfun(@ge, y, b.'), 2);
    dataRX = bin2gray(dataRX, 'pam', M).';

    % True BER
    [~, ber] = biterr(dataRX, dataTX);

    [BERtarget, ber]
    
    function r = intersection_point(a1, a2)
        r = roots([sigk(a1)^2 - sigk(a2)^2,...
            2*(a1*sigk(a2)^2 - a2*sigk(a1)^2),...
            a2^2*sigk(a1)^2 - a1^2*sigk(a2)^2 + 2*sigk(a1)^2*sigk(a2)^2*log(sigk(a2)/sigk(a1))]);

        r(r <= a1 | imag(r) ~= 0) = [];

        if isempty(r)
            warning('no solution')
            r = a1;
        end
    end
end