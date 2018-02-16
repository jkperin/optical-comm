function [x, F, exitflag, gradF, H] = saddle_free_newton(FUN, x0, options)
%% Optimize FUN using the saddle-free Newton method. The Hessian is estimated using finite difference of the gradient
% using the Gauss-Newton approximation for the Hessian

mu = options.AdaptationConstant; 
h = options.FiniteDiffStepSize;
N = options.MaxIterations;
AbsTol = options.AbsTol;
MinStep = options.MinStep;

x = x0;
[F, gradF] = FUN(x);
fprintf('0: F = %.3f\n', FUN(x))
iteration = 1;
tol = Inf;
stepSize = Inf;
prev_F = Inf;
while iteration < N && tol > AbsTol && stepSize > MinStep 
    % FUN must return value and gradient
    H = zeros(length(x0));
    parfor n = 1:length(x0)
        xh = x;
        xh(n) = xh(n) + h;
        [~, gradFh] = FUN(xh);
        H(:, n) = (gradFh - gradF)/h; % Hessian estimate
    end
    
    % Diagonalize H
    [V, D] = eig(H);
    absHinv = V*diag(abs(1./diag(D)))*V.';
    
    dx = -mu*(absHinv*gradF); % absolute of Hessian characterizes saddle-free algorithm
   
    x = x + dx.';
    [F, gradF] = FUN(x);
    
    %
    stepSize = norm(dx, 2);
    tol = abs(F - prev_F);
    prev_F = F;
    fprintf('%d: F = %.3f, ||dx||2 = %f, ||dx||inf = %f, abs tol = %G\n', iteration, F, stepSize, norm(dx, Inf), tol)
    iteration = iteration + 1;
end

if tol <= AbsTol
    fprintf('saddle_free_newton: optimization stopped because change in objective function is smaller than absolute tolerance\n')
    exitflag = 0;
end
if stepSize <= MinStep
    fprintf('saddle_free_newton: optimization stopped because step size is smaller than minimum allowed\n')
    exitflag = 1;
end
if iteration >= N
    fprintf('saddle_free_newton: maximum number of iterations exceeded\n')
    exitflag = 2;
end
