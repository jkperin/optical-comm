function D = GN_model_coeff(lamb, Df, Fiber, l)
    %% Compute Gaussian noise (GN) model coefficients
    % These coefficients assume number of spans = 1
    % Inputs:
    % - lamb: vector of wavelengths. Wavelengths are given in m
    % - Df: channel spacing in Hz
    % - Fiber: instance of class fiber correspoding to one span
    % - l (optional, default = 0): either -1, 0, or 1. l = 0 accounts for
    % the dominant contributions, while l = +- 1 accounts for corner
    % contributions
    % Output:
    % - D: 2N-1 x 2N-1 matrix representation of the nonlinear coefficients 
    % D(n1, n2, n). D(i, j) is indexed such that i = n1-n, j = n2-n.

    N = length(lamb); % Number of channels
    a = log(10)/1e4*Fiber.att(mean(lamb)); % power attenuation in (1/m)
    % Note: In Poggiolini's work, "alpha" denotes amplitude attenuation, 
    % and not power attenuation. Hence, a = 2*alpha
    b2 = Fiber.beta2(mean(lamb)); % beta2
    L = Fiber.L; % fiber length in m
       
    if not(exist('l', 'var'))
        l = 0;
    end
   
    % Compute coefficients of first quadrant of NL coefficients matrix
    D1 = zeros(N, N);
    for jj = 1:N 
        for ii = 1:N % index ii is only to allow parfor 
            if ii > jj % only compute entries below main diag since D1 is symmetric
                continue
            end
            i = ii - 1; % i = (n1-n) = (f1-f)/Df
            j = jj - 1; % j = (n2-n) = (f2-f)/Df
            fprintf('i = %d, j = %d\n', i, j)
            D1(ii, jj) = integral3(@ (x, y, z) Dfun(x, y, z, i, j, l,...
                Df, a, b2, L),...
                -Df/2, Df/2, -Df/2, Df/2, -Df/2, Df/2,... 
                'method', 'tiled', 'AbsTol', 1e-3, 'RelTol', 1e-3); % integration method 
        end
    end
    
    D1 = Fiber.gamma^2*16/27*D1; 

    % Fill out entries above main diagonal of D1 (D1 is symmetric)
    D1 = D1 + D1.' - diag(diag(D1));
    D1 = flipud(D1); % flipup corrects for indexing D1 with indeces 1,..., N as opposed to 0,..., N-1
     
    % Matrix D1 corresponds to first quadrant of matrix D 
    D = [fliplr(D1(:, 2:end)) D1];
    D = [D; flipud(D(1:end-1, :))];
end

function D = Dfun(x, y, z, i, j, l, Df, a, b2, L)
    %% GN model coefficient at particular channels
    % This function assumes rectangular spectral shape i.e., g(f) = 1, 
    % from -Df/2 to Df/2 for any channel       
    % fprod = (f1-f).*(f2-f);
    fprod = (i*Df + (x-z)).*(j*Df + (y-z));

    rho = abs((1 - exp(-a*L + 1j*4*pi^2*b2*L*fprod))./(a - 1j*4*pi^2*b2*fprod)).^2;
    % chi = 1 for Nspans = 1
    
    D = rho*(1/Df^3);
    % Note: (1/Df^3) term is necessary because spectral shape has unit area
        
    D(abs(x+y-z+l*Df) > Df/2) = 0; % zero when frequencies fall out of band
end

% function p = rho(fprod, a, b2, L)
%     %% Coefficient rho
%     p = abs((1 - exp(-a*L + 1j*4*pi^2*b2*L*fprod))./(a - 1j*4*pi^2*b2*fprod)).^2;
% end
% 
% function x = chi(fprod, Nspans, b2, L)
%     %% Coefficient chi
%     theta = 2*pi^2*fprod*b2*L;
%     x = (sin(Nspans*theta)./sin(theta)).^2;
% 
%     % Enforces limit of chi when theta -> 0
%     idx = isnan(x);
%     x(idx) = theta(idx).^2;
% end