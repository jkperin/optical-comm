function [ber, GdB] = mpam_ber_vs_Ptx(mpam, tx, fiber1, soa1, apd1, rx, sim)

sim.Ndiscard = 16; % number of symbols to be discarded from the begning and end of the sequence

%% Time and frequency
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'

dt = 1/sim.fs;
t = (0:dt:(sim.N-1)*dt).';
df = 1/(dt*sim.N);
f = (-sim.fs/2:df:sim.fs/2-df).';

sim.t = t;
sim.f = f;

% Modulator frequency response
tx.modulator.fc = tx.modulator.Fc(1); % modulator cut off frequency
tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
tx.modulator.h = @(t) [0*t(t < 0) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0))];
tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

if isfield(tx, 'RIN_shape') 
    RIN_shapedB = tx.RIN_shape(f/1e9, tx.modulator.fc/1e9);
    RIN_shapedB = RIN_shapedB*tx.RIN_variation/max(RIN_shapedB) + tx.RIN - tx.RIN_variation;
    RIN_shape = 10.^(RIN_shapedB/10 + realmin);
    RIN_shapedB = 10*log10(10.^(tx.RIN/10)*tx.RIN_bw*RIN_shape/trapz(sim.f(abs(sim.f) <= tx.RIN_bw), RIN_shape(abs(sim.f) <= tx.RIN_bw)));
    tx.RIN_shapedB = RIN_shapedB;
    RIN_shape = 10.^(RIN_shapedB/10 + realmin);
    
    persistent RIN_plot;

    if isempty(RIN_plot) || ~isvalid(RIN_plot)
        RIN_plot = figure;
        box on, hold on, grid on
    else
        figure(RIN_plot)
    end
    plot(f/1e9, RIN_shapedB);
    plot([0 tx.RIN_bw/1e9], [1 1]*10*log10(trapz(sim.f(abs(sim.f) <= tx.RIN_bw), RIN_shape(abs(sim.f) <= tx.RIN_bw))/(2*tx.RIN_bw)))
    xlabel('Frequency (GHz)')
    ylabel('RIN (dB)')
    legend('RIN(f)', 'Equivalent White Noise')
    title('RIN Power Spectrum Density')
    axis([0 tx.RIN_bw/1e9 tx.RIN-tx.RIN_variation tx.RIN+tx.RIN_variation]);
end

%% Receiver
if strcmp(rx.filter.type, 'matched')
    % Electric Lowpass Filter
    rx.elefilt = design_filter('matched', @(t) mpam.pshape(t), 1/sim.Mct);
else
    rx.elefilt = design_filter(rx.filter.type, rx.filterN, rx.filterBw/(sim.fs/2));
end

GdB = [];
if ~isempty(soa1)
    % KLSE Fourier Series Expansion (done here because depends only on filters
    % frequency response)
    % klse_fourier(rx, sim, N, Hdisp)
    [rx.U_fourier, rx.D_fourier, rx.Fmax_fourier] = klse_fourier(rx, sim, sim.Mct*(mpam.M^sim.L + 2*sim.L)); 
    
    % BER
    disp('BER with SOA')
    [ber, mpam] = soa_ber(mpam, tx, fiber1, soa1, rx, sim);
    1;
elseif ~isempty(apd1)
    [ber, mpam, apd1] = apd_ber(mpam, tx, fiber1, apd1, rx, sim);      
    
    GdB = apd1.GaindB;
    ber.est = ber.gauss;
else
    pin = apd(0, 0, Inf, rx.R, rx.Id);
    
    [ber, mpam] = apd_ber(mpam, tx, fiber1, pin, rx, sim);
    ber.est = ber.gauss;
end

if strcmp(mpam.level_spacing, 'optimized')
    persistent levels_plot;

    if isempty(levels_plot) || ~isvalid(levels_plot)
        levels_plot = figure;
        box on, hold on, grid on
    else
        figure(levels_plot)
        hold on
    end
    
    mpam.norm_levels;
    mpam.a = mpam.a*(mpam.M-1);
    mpam.b = mpam.b*(mpam.M-1);
    plot(-10, -10, 'k', 'LineWidth', 2)
    plot(-10, -10, '--k', 'LineWidth', 2)
    plot([-0.5 0.5], mpam.a*[1 1], 'k', 'LineWidth', 2)    
    plot([-0.5 0.5], mpam.b*[1 1], '--k', 'LineWidth', 2)    
    axis([-0.5 0.5 0 mpam.M])
    legend('Levels', 'Decision Thresholds')
    hold off
end
    

