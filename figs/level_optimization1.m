function [a,b] = level_optimization1(M, BERtarget, rexdB, var1, vars_ratio)
    %% Algorithm 1
    % Error probability under a single tail for a given symbol
    Pe = log2(M)*BERtarget*(M/(2*(M-1)));

    % Noises variance
    % varTherm = rx.N0*rx.elefilt.noisebw(sim.fs)/2; % variance of thermal noise
    % Shot noise variance = Agrawal 4.4.17 (4th edition)
    % calc_noise_std = @(Plevel) sqrt(varTherm + 2*apd.q*apd.Gain^2*apd.Fa*(apd.R*Plevel/apd.Gain + apd.Id)*rx.elefilt.noisebw(sim.fs)/2);
    var2 = @(Plevel) (var1*vars_ratio)*Plevel;

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
            sig = sqrt(var1 + var2(a(level)));

            [dPthresh, ~, exitflag] = fzero(@(dPthresh) qfunc(abs(dPthresh)/sig) - Pe, 0);

            if exitflag ~= 1
                warning('level_spacing_optm: threshold optimization did not converge');
            end

            b(level) = a(level) + abs(dPthresh);

            % Find next level  
            [dPlevel, ~, exitflag] = fzero(@(dPlevel) qfunc(abs(dPlevel)/sqrt(var1 + var2(b(level) + abs(dPlevel)))) - Pe, 0);    

            if exitflag ~= 1
                warning('level_spacing_optm: level optimization did not converge');     
            end

            a(level+1) = b(level) + abs(dPlevel);
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
end