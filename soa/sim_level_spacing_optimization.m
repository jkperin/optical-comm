%% Simulations for SOA with noneq_spaced and eq_spaced level spacing compared with Gaussian approximation
function sim_level_spacing_optimization
close all, clc
format compact

addpath ../f    % general functions
addpath ../mpam
addpath ../soa
addpath ../soa/f

M = [2 4];
FndB = 10:-1:3;

Colors = {'k', 'b', 'r', 'g'};
figure(1), hold on, grid on, box on
for m = 1:length(M)
    [Ptx{m}, ber{m}, Gsoa{m}] = calc_power_sensitivity(M(m), FndB);
    
    figure(1)
    plot(FndB, Ptx{m}.eq_spaced, '-', 'Color', Colors{m})
    plot(FndB, Ptx{m}.optimized, '--', 'Color', Colors{m})
    plot(FndB, Ptx{m}.optimized_gauss, ':', 'Color', Colors{m})
    1;
    
%     save partial_results_2and4PAM M FndB Ptx ber Gsoa
end
legend('Uniform Level Spacing', 'Non-Uniform Level Spacing', 'Non-Uniform Level Spacing with Gaussian Approximation', 'Location', 'NorthWest')
xlabel('Noise Figure (dB)')
ylabel('Transmitted Power (dBm)')

% save results_2and4PAM M FndB Ptx ber Gsoa

end

function [Ptx, ber, Gsoa] = calc_power_sensitivity(M, FndB)
    % Simulation parameters
    sim.Nsymb = 2^15; % Number of symbols in montecarlo simulation
    sim.Mct = 15;     % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
%     sim.L = 2;        % de Bruijin sub-sequence length (ISI symbol length)   
    sim.Me = 16; % number of used eigenvalues
    sim.BERtarget = 1e-4; 
    sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
    sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

    sim.optimize_gain = false;
    
    sim.shot = false; % include shot noise in montecarlo simulation (always included for pin and apd case)
    sim.RIN = false; % include RIN noise in montecarlo simulation
    sim.verbose = false; % show stuff

    % M-PAM
    % M, Rb, leve_spacing, pshape
    mpam = PAM(M, 100e9, 'equally-spaced', @(n) double(n >= 0 & n < sim.Mct));

    %% Time and frequency
    sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'

    dt = 1/sim.fs;
    t = (0:dt:(sim.N-1)*dt).';
    df = 1/(dt*sim.N);
    f = (-sim.fs/2:df:sim.fs/2-df).';

    sim.t = t;
    sim.f = f;

    %% Transmitter
    switch mpam.M
        case 2
            sim.optimize_gain = false;
            tx.PtxdBm = -35:2:-24;
            sim.L = 3;
        case 4
            tx.PtxdBm = -26:2:-10;
            sim.L = 3;
        case 8
            tx.PtxdBm = -22:2:-4;
            sim.L = 3;        
        case 16
            tx.PtxdBm = -18:2:0;
            sim.L = 1;
    end

    tx.lamb = 1310e-9; % wavelength
    tx.alpha = 0; % chirp parameter
    tx.RIN = -150;  % dB/Hz
    tx.rexdB = -15;  % extinction ratio in dB. Defined as Pmin/Pmax

    % Modulator frequency response
    tx.kappa = 1; % controls attenuation of I to P convertion

    %% Fiber
    b2b = fiber(); % fiber(L, att(lamb), D(lamb))

    %% Receiver
    rx.N0 = (20e-12).^2; % thermal noise psd
    rx.Id = 10e-9; % dark current
    rx.R = 1; % responsivity
    % Electric Lowpass Filter
%     rx.elefilt = design_filter('bessel', 5, mpam.Rs/(sim.fs/2));
    rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);
    % Optical Bandpass Filter
    rx.optfilt = design_filter('fbg', 0, 200e9/(sim.fs/2));

    % KLSE Fourier Series Expansion (done here because depends only on filters
    % frequency response)
    % klse_fourier(rx, sim, N, Hdisp)
    [rx.U_fourier, rx.D_fourier, rx.Fmax_fourier] = klse_fourier(rx, sim, sim.Mct*(mpam.M^sim.L + 2*sim.L)); 

    %% SOA
    % soa(GaindB, NF, lambda, maxGaindB)
    soa1 = soa(10, 9, 1310e-9, 20); 
    soa2 = soa(10, 9, 1310e-9, 20); 
    soa3 = soa(10, 9, 1310e-9, 20); 
   
    figure, hold on, grid on, box on   
    for k = 1:length(FndB)
        fprintf('----- Fn = %d dB -----\n', FndB(k))
        soa1.Fn = FndB(k);
        soa2.Fn = FndB(k);
        soa3.Fn = FndB(k);
        
        %% Equally-spaced levels
        disp('Equally-spaced levels')
        mpam.level_spacing = 'equally-spaced';
       
        if sim.optimize_gain && ~(soa1.GaindB >= 0.97*soa1.maxGaindB)
            soa1.optimize_gain(mpam, tx, b2b, rx, sim)
        else
            soa1.GaindB = soa1.maxGaindB;
        end
        Gsoa.eq_spaced(k) = soa1.GaindB;
        soa1.GaindB

        ber.eq_spaced(k) = soa_ber(mpam, tx, b2b, soa1, rx, sim);

        Ptx.eq_spaced(k) = interp1(log10(ber.eq_spaced(k).est), tx.PtxdBm, log10(sim.BERtarget), 'spline');
       
        %% Optimized level spacing
        disp('Optimized level spacing')
        sim.stats = 'accurate';
        mpam.level_spacing = 'optimized';
        
        if sim.optimize_gain && ~(soa2.GaindB >= 0.97*soa2.maxGaindB)
            soa2.optimize_gain(mpam, tx, b2b, rx, sim)
        else
            soa2.GaindB = soa2.maxGaindB;
        end
        Gsoa.optimized(k) = soa2.GaindB;
        soa2.GaindB
        
        ber.optimized(k) = soa_ber(mpam, tx, b2b, soa2, rx, sim);

        Ptx.optimized(k) = interp1(log10(ber.optimized(k).est), tx.PtxdBm, log10(sim.BERtarget), 'spline'); 
                
        %% Non-eq_spaced level spacing with Gaussian approximation
        disp('Optimized level spacing with Gaussian approximation')
        sim.stats = 'gaussian';
        mpam.level_spacing = 'optimized';
        
        if sim.optimize_gain && ~(soa3.GaindB >= 0.97*soa3.maxGaindB)
            soa3.optimize_gain(mpam, tx, b2b, rx, sim)
        else
            soa3.GaindB = soa3.maxGaindB;
        end
        Gsoa.optimized_gauss(k) = soa3.GaindB;
        soa3.GaindB
        
        ber.optimized_gauss(k) = soa_ber(mpam, tx, b2b, soa3, rx, sim);

        Ptx.optimized_gauss(k) = interp1(log10(ber.optimized_gauss(k).est), tx.PtxdBm, log10(sim.BERtarget), 'spline'); 

        plot(tx.PtxdBm, log10(ber.eq_spaced(k).est), '-b')
        plot(tx.PtxdBm, log10(ber.eq_spaced(k).count), '--ob')

        plot(tx.PtxdBm, log10(ber.optimized(k).est), '-r')
        plot(tx.PtxdBm, log10(ber.optimized(k).count), '--or')
        
        plot(tx.PtxdBm, log10(ber.optimized_gauss(k).est), '-g')
        plot(tx.PtxdBm, log10(ber.optimized_gauss(k).count), '--og')
    end

    xlabel('Received Power (dBm)')
    ylabel('log(BER)')
    legend('KLSE Fourier', 'Montecarlo', 'Location', 'SouthWest')
    axis([tx.PtxdBm(1) tx.PtxdBm(end) -8 0])
    set(gca, 'xtick', tx.PtxdBm)
    saveas(gca, sprintf('results_%dPAM.png', mpam.M))

    %% Figures
%     figure, hold on, grid on, box on
%     plot(FndB, Prx_eq_spaced)
%     plot(FndB, Prx_noneq_spaced)
%     legend('Uniform Level Spacing', 'Non-Uniform Level Spacing', 'Location', 'NorthWest')
%     xlabel('Noise Figure (dB)')
%     ylabel('Receiver Sensitivity (dBm)')
end

%% Plot Frequency response
% signal = design_filter('matched', mpam.pshape, 1/sim.Mct);
% Hsig = signal.H(sim.f/sim.fs); % signal frequency response
% figure, box on, grid on, hold on
% plot(f/1e9, abs(Hsig).^2)
% if isfield(tx, 'modulator')
%     plot(f/1e9, abs(tx.modulator.H(f)).^2)
% else
%     plot(f/1e9, ones(size(f)))
% end
% plot(f/1e9, abs(fiber.Hfiber(f, tx)).^2)
% plot(f/1e9, abs(rx.optfilt.H(f/sim.fs)).^2)
% plot(f/1e9, abs(rx.elefilt.H(f/sim.fs)).^2)
% legend('Signal', 'Modulator', 'Fiber frequency response (small-signal)', 'Optical filter', 'Receiver electric filter')
% xlabel('Frequency (GHz)')
% ylabel('|H(f)|^2')
% axis([0 rx.Fmax_fourier*sim.fs/1e9 0 3])
    