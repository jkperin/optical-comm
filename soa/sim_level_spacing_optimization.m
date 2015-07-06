%% Simulations for SOA with nonuniform and uniform level spacing compared with Gaussian approximation
function sim_level_spacing_optimization
close all, clc
format compact

addpath ../f    % general functions
addpath ../mpam/f
addpath ../soa
addpath ../soa/f

M = [2 4];
FndB = 10:-1:3;

Colors = {'k', 'b', 'r', 'g'};
figure(1), hold on, grid on, box on
for m = 1:length(M)
    [Ptx{m}, ber{m}, Gsoa{m}] = calc_power_sensitivity(M(m), FndB);
    
    figure(1)
    plot(FndB, Ptx{m}.uniform, '-', 'Color', Colors{m})
    plot(FndB, Ptx{m}.nonuniform, '--', 'Color', Colors{m})
    plot(FndB, Ptx{m}.nonuniform_gauss_approx, ':', 'Color', Colors{m})
    1;
    
    save partial_results_2and4PAM M FndB Ptx ber Gsoa
end
legend('Uniform Level Spacing', 'Non-Uniform Level Spacing', 'Non-Uniform Level Spacing with Gaussian Approximation', 'Location', 'NorthWest')
xlabel('Noise Figure (dB)')
ylabel('Transmitted Power (dBm)')

save results_2and4PAM M FndB Ptx ber Gsoa

end

function [Ptx, ber, Gsoa] = calc_power_sensitivity(M, FndB)
    % Simulation parameters
    sim.Nsymb = 2^15; % Number of symbols in montecarlo simulation
    sim.Mct = 15;     % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
%     sim.L = 2;        % de Bruijin sub-sequence length (ISI symbol length)
    sim.M = 4; % Ratio of optical filter BW and electric filter BW (must be integer)
    sim.Me = 16; % number of used eigenvalues
    sim.verbose = ~true; % show stuff
    sim.BERtarget = 1e-4; 
    sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
    sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

    sim.optimize_gain = true;
    
    sim.shot = false; % include shot noise in montecarlo simulation (always included for pin and apd case)
    sim.RIN = false; % include RIN noise in montecarlo simulation

    % M-PAM
    mpam.M = M;
    mpam.Rb = 100e9;
    mpam.Rs = mpam.Rb/log2(mpam.M);
    mpam.pshape = @(n) double(n >= 0 & n < sim.Mct); % pulse shape

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
    rx.optfilt = design_filter('fbg', 0, sim.M*rx.elefilt.fcnorm);

    % KLSE Fourier Series Expansion (done here because depends only on filters
    % frequency response)
    % klse_fourier(rx, sim, N, Hdisp)
    [rx.U_fourier, rx.D_fourier, rx.Fmax_fourier] = klse_fourier(rx, sim, sim.Mct*(mpam.M^sim.L + 2*sim.L)); 

    %% SOA
    % soa(GaindB, NF, lambda, maxGaindB)
    soa1 = soa(10, 9, 1310e-9, 20); 
    soa2 = soa(10, 9, 1310e-9, 20); 
    soa3 = soa(10, 9, 1310e-9, 20); 
    
    % uniform level spacing
    mpam_uniform = mpam; 
    mpam_uniform.level_spacing = 'uniform';
    
    mpam_nonuniform = mpam; 
    mpam_nonuniform.level_spacing = 'nonuniform';
    mpam_nonuniform.level_spacing_with_gaussian_approx = false;
    
    mpam_nonuniform_gauss = mpam; 
    mpam_nonuniform_gauss.level_spacing = 'nonuniform';
    mpam_nonuniform_gauss.level_spacing_with_gaussian_approx = true;
    
    clear mpam

    figure, hold on, grid on, box on   
    for k = 1:length(FndB)
        fprintf('----- Fn = %d dB -----\n', FndB(k))
        soa1.Fn = FndB(k);
        soa2.Fn = FndB(k);
        soa3.Fn = FndB(k);
        
        %% Uniform level spacing
        disp('Uniform level spacing')
       
        if sim.optimize_gain && ~(soa1.GaindB >= 0.97*soa1.maxGaindB)
            soa1.optimize_gain(mpam_uniform, tx, b2b, rx, sim)
        else
            soa1.GaindB = soa1.maxGaindB;
        end
        Gsoa.uniform(k) = soa1.GaindB;
        soa1.GaindB

        ber.uniform(k) = soa_ber(mpam_uniform, tx, b2b, soa1, rx, sim);

        Ptx.uniform(k) = interp1(log10(ber.uniform(k).est), tx.PtxdBm, log10(sim.BERtarget), 'spline');
       
        %% Non-uniform level spacing
        disp('Non-uniform level spacing')
        
        if sim.optimize_gain && ~(soa2.GaindB >= 0.97*soa2.maxGaindB)
            soa2.optimize_gain(mpam_nonuniform, tx, b2b, rx, sim)
        else
            soa2.GaindB = soa2.maxGaindB;
        end
        Gsoa.nonuniform(k) = soa2.GaindB;
        soa2.GaindB
        
        ber.nonuniform(k) = soa_ber(mpam_nonuniform, tx, b2b, soa2, rx, sim);

        Ptx.nonuniform(k) = interp1(log10(ber.nonuniform(k).est), tx.PtxdBm, log10(sim.BERtarget), 'spline'); 
                
        %% Non-uniform level spacing with Gaussian approximation
        disp('Non-uniform level spacing with Gaussian approximation')
        
        if sim.optimize_gain && ~(soa3.GaindB >= 0.97*soa3.maxGaindB)
            soa3.optimize_gain(mpam_nonuniform_gauss, tx, b2b, rx, sim)
        else
            soa3.GaindB = soa3.maxGaindB;
        end
        Gsoa.nonuniform_gauss_approx(k) = soa3.GaindB;
        soa3.GaindB
        
        ber.nonuniform_gauss_approx(k) = soa_ber(mpam_nonuniform_gauss, tx, b2b, soa3, rx, sim);

        Ptx.nonuniform_gauss_approx(k) = interp1(log10(ber.nonuniform_gauss_approx(k).est), tx.PtxdBm, log10(sim.BERtarget), 'spline'); 

        plot(tx.PtxdBm, log10(ber.uniform(k).est), '-b')
        plot(tx.PtxdBm, log10(ber.uniform(k).count), '--ob')

        plot(tx.PtxdBm, log10(ber.nonuniform(k).est), '-r')
        plot(tx.PtxdBm, log10(ber.nonuniform(k).count), '--or')
        
        plot(tx.PtxdBm, log10(ber.nonuniform_gauss_approx(k).est), '-g')
        plot(tx.PtxdBm, log10(ber.nonuniform_gauss_approx(k).count), '--og')
    end

    xlabel('Received Power (dBm)')
    ylabel('log(BER)')
    legend('KLSE Fourier', 'Montecarlo', 'Location', 'SouthWest')
    axis([tx.PtxdBm(1) tx.PtxdBm(end) -8 0])
    set(gca, 'xtick', tx.PtxdBm)
    saveas(gca, sprintf('results_%dPAM.png', mpam_uniform.M))

    %% Figures
%     figure, hold on, grid on, box on
%     plot(FndB, Prx_uniform)
%     plot(FndB, Prx_nonuniform)
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
    