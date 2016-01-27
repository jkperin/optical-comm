  % adeq.m    Joseph M. Kahn  9-11-11
% Time-domain linear equalizer adapted by LMS algorithm.

function [xhat, xD, eps, epsD, WTe, WTo] = LMS(x, tsamp, Ysq, Nsymb, oversamp, AdEq)

Lfilt    = AdEq.LFilt;
muvec    = AdEq.muVec;
decwtvec = AdEq.decwtvec;

% initialize digital filters (one for even sample intervals, another for odd sample intervals (needed only if M/K = 3/2).
WTe = zeros(2,4*Lfilt+2);
WTo = zeros(2,4*Lfilt+2);

% initialize filter output vectors
x1hat = zeros(1,Nsymb);
x2hat = zeros(1,Nsymb);

% initialize decision vectors
x1D = zeros(1,Nsymb);
x2D = zeros(1,Nsymb);

% initialize error vectors
% based on transmitted symbols
eps1 = zeros(1,Nsymb);
eps2 = zeros(1,Nsymb);
% based on decision symbols
eps1D = zeros(1,Nsymb);
eps2D = zeros(1,Nsymb);

% split into I and Q
x1 = x(1,:);
x2 = x(2,:);
Y1sq = Ysq(1,:);
Y2sq = Ysq(2,:);

% Run LMS adaptive algorithm. The index k runs over symbols.
for k = 1:Nsymb;
    % form the Yk
    eqindk = mod(floor(oversamp*(k-1))+Lfilt:-1:floor(oversamp*(k-1))-Lfilt,size(tsamp,2))+1;% indices for Yk
    Yk = [Y1sq(eqindk),Y2sq(eqindk)].';        % form the Yk
    % select the correct WT
    if mod(oversamp*k,1) == 0                        % this case occurs when M/K = 1, 2 for all k or for M/K = 3/2 for even k
        WT = WTe;
    elseif mod(oversamp*k,1) == 1/2;                % this case occurs only when M/K = 3/2 for odd k
        WT = WTo;
    end
    % compute and store equalizer output at time k
    xkhat = WT*Yk;          % equalizer output at time k
    %xkhat = -2*W'*Yk;          % equalizer output at time k
    x1hat(k) = xkhat(1,:);  % store equalizer output 1 at time k
    x2hat(k) = xkhat(2,:);  % store equalizer output 2 at time k

    % make and store decisions at time k
    xkD = sign(real(xkhat))+1i*sign(imag(xkhat));                           % decision at time k
    x1D(k) = xkD(1,:);                                                      % store decision on output 1 at time k
    x2D(k) = xkD(2,:);                                                      % store decision on output 2 at time k

    % Compute and store error at time k. Sign of error is such that when there is no ISI, error = noise

    % based on transmitted symbols
    eps1(k) = xkhat(1,:)-x1(k);
    eps2(k) = xkhat(2,:)-x2(k);

    % based on decision symbols
    eps1D(k) = xkhat(1,:)-xkD(1,:);
    eps2D(k) = xkhat(2,:)-xkD(2,:);

    % update WT
    % step size parameter at time k
    muk = muvec(k);
    % error at time k (based on transmitted symbols or decision symbols, depending on value of decwtvec
    epsk = decwtvec(k)*(xkhat-[x1(k);x2(k)])+~decwtvec(k)*(xkhat-xkD);
    % update WT
    WT = WT - 2*muk*([epsk(1)*Yk';epsk(2)*Yk']);
    % store WT to the correct place
    if mod(oversamp*k,1) == 0                       % this case occurs when M/K = 1, 2 for all k or for M/K = 3/2 for even k
        WTe = WT;
    elseif mod(oversamp*k,1) == 1/2;                % this case occurs only when M/K = 3/2 for odd k
        WTo = WT;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WTe = -W';
% Compile outputs 
xhat = [x1hat; x2hat];
xD = [x1D; x2D];
eps = [eps1; eps2];
epsD = [eps1D; eps2D];

end