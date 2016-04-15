function [X, Xhat] = dpll(Y, DPLL, sim)
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
if strcmpi(DPLL.CT2DT, 'bilinear')
    wn = 2/Ts*atan(DPLL.wn*Ts/2); % Compensates for frequency warping in bilinear transformation.
else
    wn = DPLL.wn;
end
Delay = DPLL.Delay; % Number of symbols in loop delay
Ntrain = DPLL.Ntrain;
Ytrain = DPLL.Ytrain;

% Generate continuous-time loop filter and convert it to discrete-time
% using bilinear transformation
% s = tf('s');    
% Fs = 2*csi*wn + wn^2/s;
% G = Fs/s;
% H = G/(1 + G); % loop filter
% [nums, dens] = tfdata(H);
% nums = cell2mat(nums);
% dens = cell2mat(dens); 
% Open-loop analog filter coefficients
nums = [0 2*csi*wn wn^2];
dens = [1 0 0]; % descending powers of s
% Closed-loop digital filter coefficients using bilinear transformation
if strcmpi(DPLL.CT2DT, 'bilinear')
    [numz, denz] = bilinear(nums, dens+nums, sim.Rs); % ascending powers of z^–1
elseif strcmpi(DPLL.CT2DT, 'impinvar')
    [numz, denz] = impinvar(nums, dens+nums, sim.Rs); % ascending powers of z^–1
else
    error('dpll/invalid method for continuous time to discrete time conversion')
end
Nnum = length(numz);
Nden = length(denz);
window_num = 0:-1:-Nnum+1;
window_den = 0:-1:-Nden+2;

N = length(Y);
X = zeros(2, N);
phiLO = zeros(2, N+1);
phiN = zeros(2, N+1);
Xhat = zeros(2, N);
for k = max([Nnum, Nden, Delay]+1):N % runs over symbols
    % Correct phase
    X(:, k) = Y(:, k).*exp(-1j*phiLO(:, k-Delay));
    
    % Phase estimation
    if strcmpi(DPLL.phaseEstimation, '4th power')
        phiN_tilde = 1/M*angle(-Y(:, k).^4);        
        % Unwrap
        p  = floor((phiN(:,k-1)-phiN_tilde+pi)/(2*pi)); % unwrap
        phiN(:, k) = phiN_tilde + p*(2*pi);        
    elseif strcmpi(DPLL.phaseEstimation, 'DD')
        if k < Ntrain % training sequence
            Xhat(:, k) = Ytrain(:, k);
        else % switch to decision directed
            Xhat(:, k) = detect(X(:, k)); 
        end
        
        % Remove signal component from phase
        phiN_tilde = angle(Y(:, k)) - angle(Xhat(:, k));
        p  = floor((phiN(:,k-1)-phiN_tilde+pi)/(2*pi)); % unwrap
        phiN(:, k+1) = phiN_tilde + p*(2*pi);
        % phiN is the unwrapped sum of phase noise and AWGN noise
        % phase component
    else
        error('dpll/invalid option for phase estimation')
    end
    
    % Loop filter updates LO phase
    phiLO(:, k+1) = phiN(:, k+1 + window_num)*numz.'...
        - phiLO(:, k + window_den)*denz(2:end).';
end
% 
phiN = unwrap(phiN, [], 2);
phiLO = unwrap(phiLO, [], 2);

% Plots
if isfield(sim, 'Plots') && sim.Plots.isKey('DPLL phase error') && sim.Plots('DPLL phase error')
    figure(200), clf
    subplot(211), box on, hold on
    plot(phiN(1, :), 'b')
    plot(phiLO(1, :), '--r')
    a = axis;
    plot([Ntrain Ntrain], a(3:4), '--k')
    legend('\phi_N = phase without signal component', '\phi_{LO} = local oscillator phase',...
        'End of training', 'Location', 'Best')
    hold off
    
    subplot(212), box on, hold on
    plot(phiN(1, :)-phiLO(1, :), 'b')
    a = axis;
    plot([Ntrain Ntrain], a(3:4), '--k')
    legend('\phi_N-\phi_{LO}', 'End of training', 'Location', 'Best')
    title(['var(\phi_s)/var(\phi_s-\phi_{LO}) = ' num2str(var(phiN(1, :))/var(phiN(1, :)-phiLO(1, :)))])
    hold off
    drawnow
end

function x = detect(y)
    x = qammod(qamdemod(y, M, 0, 'Gray'), M, 0, 'Gray');
end
end

