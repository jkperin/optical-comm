function margin_vs_gain()
clc, close all

addpath ../mpam
addpath ../f
addpath f

% Simulation parameters
sim.Nsymb = 2^14; % Number of symbols in montecarlo simulation
sim.Mct = 9;    % Oversampling ratio to simulate continuous time (must be odd so that sampling is done  right, and FIR filters have interger grpdelay)  
sim.L = 2;        % de Bruijin sub-sequence length (ISI symbol length)
sim.BERtarget = 1.8e-4; 
sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation

%
sim.shot = true; % include shot noise. Only included in montecarlo simulation (except for APD)
sim.RIN = true; % include RIN noise. Only included in montecarlo simulation
sim.verbose = false; % show stuff

% M-PAM
mpam = PAM(4, 100e9, 'equally-spaced', @(n) double(n >= 0 & n < sim.Mct));

%% Time and frequency
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'

dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

sim.t = t;
sim.f = f;

%% Transmitter
tx.lamb = 1310e-9; % wavelength
tx.alpha = 0; % chirp parameter
tx.RIN = -150;  % dB/Hz
tx.rexdB = -10;  % extinction ratio in dB. Defined as Pmin/Pmax

tx.PtxdBm = -22:1:-8; % range for equally spaced levels

% Modulator frequency response
% tx.modulator.fc = 30e9; % modulator cut off frequency
% tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
% tx.modulator.h = @(t) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0));
% tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

%% b2b
b2b = fiber();

%% Receiver
rx.N0 = (30e-12).^2; % thermal noise psd
% Electric Lowpass Filter
% rx.elefilt = design_filter('bessel', 5, 1.25*mpam.Rs/(sim.fs/2));
rx.elefilt = design_filter('matched', mpam.pshape, 1/sim.Mct);
rx.Id = 10e-9;
rx.R = 1;

%% Equalization
rx.eq.type = 'None';
% rx.eq.ros = 2;
rx.eq.Ntaps = 31;
% rx.eq.Ntrain = 2e3;
% rx.eq.mu = 1e-2;

%% PIN
pin = apd(0, 0, Inf, rx.R, rx.Id);
ber_pin = apd_ber(mpam, tx, b2b, pin, rx, sim);

% Variables to iterate
Gains = 1:20;
ka = [0.1 0.25 0.5];
BW = 1e9*[20 100; 20 200; 20 300;...
    30 100; 30 200; 30 300]; % (10:2.5:50)*1e9;

MargindB = zeros(length(ka), length(BW), length(Gains));
Gopt_margin = zeros(length(ka), length(BW));
OptMargindB = zeros(length(ka), length(BW));

legs = {};
for n = 1:length(ka)
    figure(1000+n), hold on, box on, grid on
    for m = 1:size(BW, 1)
        try
            [BER{n, m}, MargindB(n, m, :), Gopt_margin(n, m), OptMargindB(n, m)]...
                = iterate(Gains, mpam, ka(n), BW(m, :), sim); 
            
            figure(1000+n), hold on
            margin = MargindB(n, m, :);
            margin = margin(:);
            hlinem(m) = plot(Gains, margin);
            plot(Gopt_margin(n, m), OptMargindB(n, m), 'o', 'Color', get(hlinem(m), 'Color'));
            legs = [legs sprintf('ka = %.2f, BW0 = %.2f, GainBW = %.2f', ka(n), BW(m, 1)/1e9, BW(m, 2)/1e9)];
        catch e
            warning(e.message)
            MargindB(n, m, :) = NaN;
            Gopt_margin(n, m) = NaN;
            OptMargindB(n, m) = NaN;
        end       
    end
    figure(1000+n)
    xlabel('APD Gain (Linear Units)')
    ylabel('Margin Improvement (dB)')
    legend(hlinem, legs, 'Location', 'SouthEast')
    drawnow
end
 
% figure, hold on, box on
% legs = {};
% for n = 1:length(ka)
%     plot(BW/1e9, OptMargindB(n, :), '-');
%     legs = [legs sprintf('ka = %.2f', ka(n))];
% end
% xlabel('Bandwidth (GHz)')
% ylabel('Optimal Margin (dB)')
% legend(legs)
% 
% 
% figure, box on, hold on
% legs = {};
% for n = 1:length(ka)
%     plot(BW/1e9, Gopt_margin(n, :), '-');
%     legs = [legs sprintf('ka = %.2f', ka(n))];
% end
% xlabel('Bandwidth (GHz)')
% ylabel('Optimal Gain (dB)')
% legend(legs)

% Save results
filename = sprintf('margin_vs_gain_BW_results_%d-PAM_%s', mpam.M, mpam.level_spacing);
sim = rmfield(sim, 't');
sim = rmfield(sim, 'f');
save(filename)

function [BER, MargindB, Gopt_margin, OptMargindB] = iterate(Gains, mpam, ka, BW, sim)
    % BER [1 x length(Gains)] struct = BER struct for each gain in Gains
    % MargindB [1 x length(Gains)] double = margin in dB for each gain in Gains
    % Gopt_margin = optimal gain (gain that achieves optimal margin)
    % OptMargindB = optimal margin in dB
    disp('--- Iterate ---')                  
    
    %% APD 
    % (GaindB, ka, BW, R, Id) 
    apdG = apd(10, ka, BW, rx.R, rx.Id);
    
    PtxdBm_pin_BERtarget = interp1(log10(ber_pin.gauss), tx.PtxdBm, log10(sim.BERtarget));
    
    % APD
    PtxdBm_BERtarget = zeros(size(Gains));
    figure, hold on, box on
    leg = {};
    
    hline(1) = plot(tx.PtxdBm, log10(ber_pin.gauss), '-b');
    plot(tx.PtxdBm, log10(ber_pin.awgn), '--', 'Color', get(hline(1), 'Color'))
    plot(tx.PtxdBm, log10(ber_pin.count), '-o', 'Color', get(hline(1), 'Color'))
    leg = [leg 'PIN'];
    for k= 1:length(Gains)
        fprintf('- %d-PAM, %s, ka = %.2f, BW0 = %.2f, GainBW = %.2f, Gain = %.2f\n', mpam.M, mpam.level_spacing, ka, BW(1)/1e9, BW(2)/1e9, Gains(k))
        
        apdG.Gain = Gains(k);

        % BER
        ber_apd = apd_ber(mpam, tx, b2b, apdG, rx, sim);
        BER(k) = ber_apd;

        % Calculate power at the target BER
        PtxdBm_BERtarget(k) = interp1(log10(ber_apd.gauss), tx.PtxdBm, log10(sim.BERtarget));
        
        % Figures
        hline(k+1) = plot(tx.PtxdBm, log10(ber_apd.gauss), '-');
        plot(tx.PtxdBm, log10(ber_apd.awgn), '--', 'Color', get(hline(k+1), 'Color'))
        plot(tx.PtxdBm, log10(ber_apd.count), '-o', 'Color', get(hline(k+1), 'Color'))
        leg = [leg sprintf('Gain = %.2f', apdG.Gain)];
    end   
    xlabel('Received Power (dBm)')
    ylabel('log(BER)') 
    axis([tx.PtxdBm(1) tx.PtxdBm(end) -8 0])
    set(gca, 'xtick', tx.PtxdBm)
    title(sprintf('%d-PAM, %s, ka = %.2f, BW = %.2f', mpam.M, mpam.level_spacing, ka, BW/1e9))
           
    % Calculate margin   
    MargindB = PtxdBm_pin_BERtarget - PtxdBm_BERtarget;

    % Find gain for optimal margin
    sim.OptimizeGain = true;
    
    % Optimal gain
    [ber_apd_opt, ~, apdG] = apd_ber(mpam, tx, b2b, apdG, rx, sim);
    Gopt_margin = apdG.Gain;
    
    PtxdBm_BERtarget_opt = interp1(log10(ber_apd_opt.gauss), tx.PtxdBm, log10(sim.BERtarget));
    
    OptMargindB = PtxdBm_pin_BERtarget - PtxdBm_BERtarget_opt;  
    
    % Plot optimal gain
    hline(end+1) = plot(tx.PtxdBm, log10(ber_apd_opt.gauss), '-k');
    plot(tx.PtxdBm, log10(ber_apd_opt.awgn), '--', 'Color', get(hline(end), 'Color'))
    plot(tx.PtxdBm, log10(ber_apd_opt.count), '-o', 'Color', get(hline(end), 'Color'))
    leg = [leg sprintf('Optimal Gain = %.2f', apdG.Gain)];
    legend(hline, leg)  
    drawnow
end
end
