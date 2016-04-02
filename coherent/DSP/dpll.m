function X = dpll(Y, DPLL, sim)
%% Digital phase locked loop
% Loop filter is assumed to be a second-order filter defined by parameters
% csi (damping) and wn (relaxation frequency). The loop filter frequency
% response is calculated in continuous time and converted to discrete time
% using biliniear transformation.
% Inputs:
% - Y : input signal in two polarizations
% - DPLL : struct containing DPLL parameters DPLL.{csi, wn, Ntrain, Ytrain}
% - sim : simulation parameters sim.{M : modulation order, Rs : symbol rate}

M = sim.M;
Ts = 1/sim.Rs;
csi = DPLL.csi;
wn = 2/Ts*atan(DPLL.wn*Ts/2); % Compensates for frequency warping in bilinear transformation.
Ntrain = DPLL.Ntrain;
Ytrain = DPLL.Ytrain;

% Generate continuous-time loop filter and convert it to discrete-time
% using bilinear transformation
% s = tf('s');    
% Fs = 2*csi*wn + wn^2/s;
% G = Fs/s;
% H = G/(1 + G); % loop filter: relates 
% [nums, dens] = tfdata(H);
% nums = cell2mat(nums);
% dens = cell2mat(dens); 
nums = [2*csi*wn wn^2];
dens = [1 0 0]; % descending powers of s
[numz, denz] = bilinear(nums, dens, sim.Rs); % ascending powers of z^–1
Nnum = length(numz);
Nden = length(denz);

N = length(Y);
X = zeros(2, N);
phiLO = zeros(2, N+1);
phiS = zeros(2, N);
for k = max(Nnum, Nden):N % runs over symbols
    % Correct phase
    X(:, k) = Y(:, k).*exp(-1j*phiLO(:, k));
    
    % Phase estimation
    if strcmpi(DPLL.phaseEstimation, '4th power')
        error('dpll/4th power not implemented yet')
%         phiS_tilde = angle(Y(:, k).^4);        
%         p  = floor((phiS(:,k-1)-phiS_tilde+pi)/(2*pi)); % unwrap
%         phiS(:, k) = phiS_tilde + p*(2*pi);        
    elseif strcmpi(DPLL.phaseEstimation, 'DD')
        if Ntrain < k % training sequence
            Xhatk = Ytrain(:, k);
        else % switch to decision directed
            Xhatk = detect(X(:, k)); 
        end
        
        phiS_tilde = angle(X(:, k)) - angle(Xhatk);
        p  = floor((phiS(:,k-1)-phiS_tilde+pi)/(2*pi)); % unwrap
        phiS(:, k) = phiS_tilde + p*(2*pi);
    else
        error('dpll/invalid option for phase estimation')
    end
    
    % Loop filter updates LO phase
    phiLO(:, k+1) = phiS(:, k:-1:k-Nnum+1)*numz.' - phiLO(:, k:-1:k-Nden+2)*denz(2:end).';
end

% figure, hold on
% plot(unwrap(angle(Y(1, :)) - angle(Ytrain(1, :))))
% plot(phiS(1, :), 'r')
% plot(phiLO(1, :), 'k')
% figure(h)

function x = detect(y)
    x = qammod(qamdemod(y, M, 0, 'Gray'), M, 0, 'Gray');
end
end

