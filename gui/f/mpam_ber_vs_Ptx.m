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
tx.kappa = 1; % controls attenuation of I to P convertion
tx.modulator.fc = tx.modulator.Fc(1); % modulator cut off frequency
tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
tx.modulator.h = @(t) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0));
tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

if isfield(tx, 'RIN_shape') 
    RIN_shapedB = tx.RIN_shape(f/1e9, tx.modulator.fc/1e9);
    RIN_shapedB = RIN_shapedB*tx.RIN_variation/max(RIN_shapedB) + tx.RIN - tx.RIN_variation;
    RIN_shape = 10.^(RIN_shapedB/10 + realmin);
    RIN_shapedB = 10*log10(10.^(tx.RIN/10)*tx.RIN_bw*RIN_shape/trapz(sim.f(abs(sim.f) <= tx.RIN_bw), RIN_shape(abs(sim.f) <= tx.RIN_bw)));
    tx.RIN_shapedB = RIN_shapedB;
    RIN_shape = 10.^(RIN_shapedB/10 + realmin);
    
    figure(5), hold on, grid on
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
    ber = soa_ber(mpam, tx, fiber1, soa1, rx, sim);
    1;
elseif ~isempty(apd1)
    if sim.OptimizeGain
        apd1.optimize_gain(mpam, tx, fiber1, rx, sim);
        GdB = apd1.GaindB;
    end
        
    ber = apd_ber(mpam, tx, fiber1, apd1, rx, sim);      
    ber.est = ber.gauss;
else
    pin = apd(0, 0, Inf, rx.R, rx.Id);
    
    ber = apd_ber(mpam, tx, fiber1, pin, rx, sim);
    ber.est = ber.gauss;
end

