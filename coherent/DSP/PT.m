% PhaseTracker     Milad Sharif  4-9-13
% Decision Directed Phase Tracker
% Kil Nam Oh, et al "Modified constant modulus algorithm: blind
% equalizaiton and carrier phase recovery algorithm", 

function [xD, z, eps, theta] = PT(x,xequ, Nsymb,Pt)


mu       = Pt.mu;
decwtvec = Pt.decwtvec;

% Initialize theta
theta =  zeros(2,Nsymb);

% equalized output with phase error correction
z1 = zeros(1,Nsymb);
z2 = zeros(1,Nsymb);

% initialize error vectors based on transmitted symbols
eps1 = zeros(1,Nsymb);
eps2 = zeros(1,Nsymb);

% Initialize filter output vectors
x1D = zeros(1,Nsymb);
x2D = zeros(1,Nsymb);

x1equ = xequ(1,:);
x2equ = xequ(2,:);

% Run LMS adaptive algorithm. The index k runs over symbols.
for k = 1:Nsymb;

    z1(k) = x1equ(k)*exp(-1j*theta(1,k));
    z2(k) = x2equ(k)*exp(-1j*theta(2,k));

    x1D(k) = decwtvec(k)*x(1,k) + ~decwtvec(k)*detect(z1(k));
    x2D(k) = decwtvec(k)*x(2,k) + ~decwtvec(k)*detect(z2(k));
    
    % Compute and store error at time k. 
    eps1(k) = z1(k)-x1D(k);
    eps2(k) = z2(k)-x2D(k);


    % Update theta
    if k < Nsymb
        theta(:,k+1) = theta(:,k) - mu*[imag(z1(k)*conj(eps1(k)));imag(z2(k)*conj(eps2(k)))];
    end
end
% Compile outputs 
xD = [x1D; x2D];
z = [z1; z2];
eps = [eps1; eps2];
end


function x = detect(y)
    x = sign(real(y))+1i*sign(imag(y));
end
