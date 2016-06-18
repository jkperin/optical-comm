function [ber, mpam] = preamplified_sys_ber(mpam, tx, Fibers, Amp, rx, sim)
%% Calculate BER of pre-amplified IM-DD system
% Inputs:
% - mpam: PAM class
% - tx: struct with transmitter paramters
% - Fibers: array of Fiber classes containing the fibers used in
% transmittion i.e., Fibers = [SMF DCF] or Fibers = [SMF];
% - Amp: pre-amplifier using SOA class
% - rx: struct with receiver parameters
% - sim: struct with simulation parameters

%% Polarizer is always disregarded so number of noise polarizations is 2
% %% Level Spacing Optimization
% if strcmp(mpam.level_spacing, 'optimized')
%     % is sim.stats is set to Gaussian, then use Gaussian approximation,
%     % otherwise uses accurate statistics
%     if isfield(sim, 'stats') && strcmp(sim.stats, 'gaussian')        
%         % Optimize levels using Gaussian approximation
%         mpam = mpam.optimize_level_spacing_gauss_approx(sim.BERtarget, tx.rexdB, noise_std, sim.verbose);     
%     else
%         % Optimize levels using accurate noise statisitics
%         [a, b] = level_spacing_optm(mpam, tx, soa, rx, sim);
%         mpam = mpam.set_levels(a, b);
%     end
% end   

% Calculate total (residual) dispersion
Dtotal = 0;
for k = 1:length(Fibers)
    fiberk = Fibers(k);
    
    Dtotal = Dtotal + fiberk.D(tx.lamb)*fiberk.L;
end

fprintf('Total dispersion at %.2f nm: %.3f ps/nm\n', tx.lamb*1e9, 1e3*Dtotal)

%% Calculations BER
% Transmitted power
Ptx = dBm2Watt(tx.PtxdBm);
ber.awgn = zeros(size(Ptx));

ber.count = zeros(size(Ptx));
for k = 1:length(Ptx)
    tx.Ptx = Ptx(k);
         
    % Montecarlo simulation
    [ber.count(k), ber.awgn(k)] = ber_preamp_sys_montecarlo(mpam, tx, Fibers, Amp, rx, sim);
    
    if sim.stopSimWhenBERreaches0 && ber.count(k) == 0
        break;
    end 
end

%% Plots
if isfield(sim, 'Plots') && sim.Plots.isKey('BER') && sim.Plots('BER') && length(ber.count) > 1
    figure(1), hold on, box on
    hline = plot(tx.PtxdBm, log10(ber.count), '-o', 'LineWidth', 2);
    plot(tx.PtxdBm, log10(ber.awgn), '-', 'Color', get(hline, 'Color'), 'LineWidth', 2)
    legend('Counted', 'AWGN approximation')
    set(gca, 'FontSize', 12)
    axis([tx.PtxdBm(1) tx.PtxdBm(end) -8 0])
end
    