function [Nsum, Nmult] = calcDSPOperations(Rx, sim)

Natan = 4 + 1; % number of real operations to calculate atan (The simplified single lookup table algorithm with nearest-neighbor linear interpolation)
Ncos = 4; % number of real operations to calculatae cos
Nsin = 4; % number of real operations to calculatae sin

complexMult = 6; % 1 complex multiplication = 6 real operations
complexAdd = 2; % 1 complex addition = 2 real operations

Npol = 2; % number of polarizations

ros = sim.ros;
[~, q] = rat(ros);

%% Equalization
eq.Nsum = 0;
eq.Nmult = 0;

if strcmpi(Rx.AdEq.structure, '2 filters')
    Nfilters = 2;
    [eq.Nsum, eq.Nmult] = countFIR(Nfilters, eq.Nsum, eq.Nmult, Rx.AdEq.Ntaps);
    [eq.Nsum, eq.Nmult] = countComplexMult(2, eq.Nsum, eq.Nmult);
    [eq.Nsum, eq.Nmult] = countComplexAdd(2, eq.Nsum, eq.Nmult);
    
elseif strcmpi(Rx.AdEq.structure, '4 filters')
    Nfilters = 4;
    [eq.Nsum, eq.Nmult] = countFIR(Nfilters, eq.Nsum, eq.Nmult, Rx.AdEq.Ntaps);
    [eq.Nsum, eq.Nmult] = countComplexAdd(2, eq.Nsum, eq.Nmult);
      
end

%% Carrier phase recovery
cpr.Nsum = 0;
cpr.Nmult = 0;
if strcmpi(Rx.CPR.type, 'DPLL')
    % Correct phase
    [cpr.Nsum, cpr.Nmult] = countComplexExp(Npol, cpr.Nsum, cpr.Nmult, 'Lookup table'); %  Y(:, k).*exp(-1j*phiLO(:, k));
    [cpr.Nsum, cpr.Nmult] = countComplexMult(Npol, cpr.Nsum, cpr.Nmult);
    
    % Phase estimation
    [cpr.Nsum, cpr.Nmult] = countComplexMult(2*Npol, cpr.Nsum, cpr.Nmult); % angle(Y(:, k)) - angle(Xhatk);
    [cpr.Nsum, cpr.Nmult] = countArctan(2*Npol, cpr.Nsum, cpr.Nmult, 'Lookup table'); 
    cpr.Nsum = cpr.Nsum + 1;
    
    [cpr.Nsum, cpr.Nmult] = countUnwrap(cpr.Nsum, cpr.Nmult); % unwrap for pol x 
    [cpr.Nsum, cpr.Nmult] = countUnwrap(cpr.Nsum, cpr.Nmult); % unwrap for pol y
    
    % Loop filter
    % Second-order DPLL closed loop response has order 5 in denominator and
    % in numerator
    Ntaps_den = 5;
    Ntaps_num = 5;
    [cpr.Nsum, cpr.Nmult] = countFIR(Npol, cpr.Nsum, cpr.Nmult, Ntaps_num+Ntaps_den-1);
    
elseif strcmpi(Rx.CPR.type, 'feedforward')
    % Correct phase
    [cpr.Nsum, cpr.Nmult] = countComplexExp(Npol, cpr.Nsum, cpr.Nmult, 'Lookup table'); %  Y(:, k).*exp(-1j*phiLO(:, k));
    [cpr.Nsum, cpr.Nmult] = countComplexMult(Npol, cpr.Nsum, cpr.Nmult);
    
    % Phase estimation
    [cpr.Nsum, cpr.Nmult] = countComplexMult(2*Npol, cpr.Nsum, cpr.Nmult); % angle(Y(:, k)) - angle(Xhatk);
    [cpr.Nsum, cpr.Nmult] = countArctan(2*Npol, cpr.Nsum, cpr.Nmult, 'Lookup table'); 
    cpr.Nsum = cpr.Nsum + 1;

    [cpr.Nsum, cpr.Nmult] = countUnwrap(cpr.Nsum, cpr.Nmult); % unwrap for pol x 
    [cpr.Nsum, cpr.Nmult] = countUnwrap(cpr.Nsum, cpr.Nmult); % unwrap for pol y
    
    % Filter
    [cpr.Nsum, cpr.Nmult] = countFIR(Npol, cpr.Nsum, cpr.Nmult, floor(0.5*(Rx.CPR.Ntaps-1))); % theta_tilde(:,k+1) = psi(:,k:-1:k-d+1) * Wsd;
    [cpr.Nsum, cpr.Nmult] = countFIR(Npol, cpr.Nsum, cpr.Nmult, Rx.CPR.Ntaps); % thetahat(k-Delta, :) = psi(:,k:-1:k-L+1) * conj(Whd);
            
    % Correct phase
    [cpr.Nsum, cpr.Nmult] = countComplexExp(Npol, cpr.Nsum, cpr.Nmult, 'Lookup table'); % x1hat(k-Delta) = y1(k-Delta)*exp(-1i*thetahat(k-Delta, 1));
    [cpr.Nsum, cpr.Nmult] = countComplexMult(Npol, cpr.Nsum, cpr.Nmult);
end

%% Phase tracking
pt.Nsum = 0;
pt.Nmult = 0;
if strcmpi(sim.ModFormat, 'QAM')
    % Phase tracker
    [pt.Nsum, pt.Nmult] = countComplexExp(Npol, pt.Nsum, pt.Nmult, 'Lookup table'); % yx_hat(n) = X(1, n)*exp(-1j*phix(n));
    [pt.Nsum, pt.Nmult] = countComplexMult(Npol, pt.Nsum, pt.Nmult);

    % Calc error
    [pt.Nsum, pt.Nmult] = countComplexAdd(Npol, pt.Nsum, pt.Nmult); % ephix = yx_hat(n) - qammod(qamdemod(yx_hat(n), M, 0, 'Gray'), M, 0, 'Gray');

    % Update coefficients
     [pt.Nsum, pt.Nmult] = countComplexMult(Npol, pt.Nsum, pt.Nmult); % phix(n+1) = phix(n) - mu*imag(yx_hat(n)*conj(ephix));
     pt.Nmult = pt.Nmult + 2*1;
     pt.Nsum = pt.Nsum + 2*1;
end


%% 
Nsum = [eq.Nsum cpr.Nsum+pt.Nsum];
Nmult = [eq.Nmult cpr.Nmult+pt.Nmult];
labels = {'Equalization', 'Phase Recovery'};
pie(Nsum+Nmult, labels)
title('Number of real fixed-point operations')

function [Nsum, Nmult] = countComplexMult(N, Nsum, Nmult)
    %% Count N complex multiplications
    Nsum = Nsum + N*2;
    Nmult = Nmult + N*4;
end

function [Nsum, Nmult] = countComplexAdd(N, Nsum, Nmult)
    %% Count N complex additions
    Nsum = Nsum + N*2;
end

function [Nsum, Nmult] = countFIR(N, Nsum, Nmult, Ntaps)
    %% Count operations in N Ntaps-FIR filters 
    Nsum = Nsum + N*(Ntaps - 1);
    Nmult = Nmult + N*(Ntaps);
end

function [Nsum, Nmult] = countUnwrap(Nsum, Nmult)
    %% Phase unwrapping operation
%     p  = floor((phiS(:,k-1)-phiS_tilde+pi)/(2*pi)); % unwrap
%     phiS(:, k) = phiS_tilde + p*(2*pi);   
    Nsum = Nsum + 4;
    Nmult = Nmult + 2;
end

function [Nsum, Nmult] = countComplexExp(N, Nsum, Nmult, algorithm)
    %% Number of real fixed-point operations to calculate exp(1j*x) = cos(x) + 1j*sin(x)
    % Reference: http://www.mathworks.com/help/fixedpoint/examples/calculate-fixed-point-sine-and-cosine.html
    % Shift operations and memory lookup operations are ignored.
    switch (algorithm)
        case 'Cordic'
            Nit = 5; % Number of iterations in cordic 
            Nsum = Nsum + N*2*3*Nit;
            % Note: 2*(+ 2 shift per iteration + 1 table lookup)
        case 'Lookup table'
            Nsum = Nsum + N*2*2;
            Nmult = Nmult + N*2*2;
            % Note: 2*(+ 2 table lookup)
    end
end
function [Nsum, Nmult] = countArctan(N, Nsum, Nmult, algorithm)
    %% Number of real fixed-point operations to calculate arctan
    % Reference: http://www.mathworks.com/help/fixedpoint/examples/calculate-fixed-point-arctangent.html
    % Shift operations and memory lookup operations are ignored.
    switch (algorithm)
        case 'Cordic'
            Nit = 5; % Number of iterations in cordic 
            Nsum = Nsum + N*3*Nit;
            % Note: + 2 shift per iteration + 1 table lookup
        case 'Chebyshev polynomial'
            Norder = 5; % Order of Chebyshev polynomial used in approximation
            Nsum = Nsum + N*(Norder - 1)/2;
            Nmult = Nmult + N*(Norder + 2);
            % Note: + 2 shift per iteration
        case 'Lookup table'
            Nsum = Nsum + 2*N;
            Nmult = Nmult + 2*N;
            % Note: + 2 table lookup
    end
end

end