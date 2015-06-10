%% Calculate BER of amplified IM-DD system. 
% BER is calculated via montecarlo simulation and through saddlepoint approximation
% to calculate tail probabilities given the moment generating funcion.

function [ber, mpam] = soa_ber(mpam, tx, soa, rx, sim)

addpath f

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
   [mpam.a, mpam.b] = level_spacing_optm(mpam, tx, soa, rx, sim);
else
    error('Invalide Option!')
end

% KLSE Fourier Series Expansion 
[U_fourier, D_fourier, Fmax_fourier] = klse_fourier(rx, sim, sim.Mct*(mpam.M^sim.L + 2*sim.L)); 

if sim.verbose
    figure, hold on
    plot(sim.f/1e9, abs(rx.elefilt.H(sim.f/sim.fs)).^2)
    plot(sim.f/1e9, abs(rx.optfilt.H(sim.f/sim.fs)).^2)
    legend('electrical', 'optical')
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
    ber.count(k) = ber_soa_montecarlo(mpam, tx, soa, rx, sim);
    
    % approx ber
    [ber.est(k), ber.gauss(k)] = ber_soa_klse_fourier(U_fourier, D_fourier, Fmax_fourier, mpam, tx, soa, rx, sim);
end

if sim.verbose
    figure, hold on
    plot(tx.PtxdBm, log10(ber.count), '-o')
    plot(tx.PtxdBm, log10(ber.est))
    plot(tx.PtxdBm, log10(ber.gauss))
    % plot(tx.PtxdBm, log10(ber.est_pdf))
    legend('Counted', 'KLSE Fourier & Saddlepoint Approx', 'Gaussian Approximation',...
        'Location', 'SouthWest')
    axis([tx.PtxdBm(1) tx.PtxdBm(end) -10 0])
end
    