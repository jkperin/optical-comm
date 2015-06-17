%% Calculate BER of unamplified IM-DD system with APD detector 
% BER is calculated via montecarlo simulation and through saddlepoint approximation
% to calculate tail probabilities given the moment generating funcion.

function [ber, mpam] = apd_ber(mpam, tx, apd, rx, sim)

dBm2Watt = @(x) 1e-3*10.^(x/10);

%% Time and frequency
dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

sim.t = t;
sim.f = f;

% Generate unipolar PAM signal
if strcmp(mpam.level_spacing, 'uniform')
    mpam.a = (0:2:2*(mpam.M-1)).';
    mpam.b = (1:2:(2*(mpam.M-1)-1)).';
elseif strcmp(mpam.level_spacing, 'nonuniform')
   [mpam.a, mpam.b] = level_spacing_optm_gauss_approx(mpam, tx, apd, rx, sim);
else
    error('Invalide Option!')
end

if sim.verbose
    figure, hold on
    plot(sim.f/1e9, abs(rx.elefilt.H(sim.f/sim.fs)).^2)
    xlabel('Frequency (GHz)')
    ylabel('|H(f)|^2')
    grid on
end

% Transmitted power
Ptx = dBm2Watt(tx.PtxdBm);

ber.count = zeros(size(Ptx));
ber.est = zeros(size(Ptx));
ber.gauss = zeros(size(Ptx));
for k = 1:length(Ptx)
    tx.Ptx = Ptx(k);
         
    % Montecarlo simulation
    ber.count(k) = ber_apd_montecarlo(mpam, tx, apd, rx, sim);
    
    % approx ber
    [~, ber.gauss(k)] = ber_apd_doubly_stochastic(mpam, tx, apd, rx, sim);
end

if sim.verbose
    figure, hold on
    plot(tx.PtxdBm, log10(ber.count), '-o')
%     plot(tx.PtxdBm, log10(ber.est))
    plot(tx.PtxdBm, log10(ber.gauss))
    % plot(tx.PtxdBm, log10(ber.est_pdf))
    legend('Counted', 'Gaussian Approximation',...
        'Location', 'SouthWest')
    axis([tx.PtxdBm(1) tx.PtxdBm(end) -10 0])
end
    