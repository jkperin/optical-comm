%% Calculate BER of amplified IM-DD system. 
% BER is calculated via montecarlo simulation and through saddlepoint approximation
% to calculate tail probabilities given the moment generating funcion.

function [ber, mpam] = soa_ber(mpam, tx, soa, rx, sim)

dBm2Watt = @(x) 1e-3*10.^(x/10);

% Generate unipolar PAM signal
if strcmp(mpam.level_spacing, 'uniform')
    mpam.a = (0:2:2*(mpam.M-1)).';
    mpam.b = (1:2:(2*(mpam.M-1)-1)).';
elseif strcmp(mpam.level_spacing, 'nonuniform')
   [mpam.a, mpam.b] = level_spacing_optm(mpam, tx, soa, rx, sim);
else
    error('Invalide Option!')
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
    [ber.est(k), ber.gauss(k)] = ber_soa_klse_fourier(rx.U_fourier, rx.D_fourier, rx.Fmax_fourier, mpam, tx, soa, rx, sim);
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
    