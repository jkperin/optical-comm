%% Power penalty vs modulator cutoff frequency
function pp = ofdm_pp_vs_mod_cutoff(ofdm, tx, fiber, rx, sim)

% Transmitter filter (ZOH + some smoothing filter)
tx.filter = design_filter('bessel', 5, 1/(ofdm.Ms*sim.Mct));

% Convolve with ZOH
bzoh = ones(1, sim.Mct)/sim.Mct;
tx.filter.num = conv(tx.filter.num, bzoh);
tx.filter.grpdelay = grpdelay(tx.filter.num, tx.filter.den, 1);
tx.filter.H = @(f) freqz(tx.filter.num, tx.filter.den, 2*pi*f).*exp(1j*2*pi*f*tx.filter.grpdelay);
tx.filter.noisebw = @(fs) noisebw(tx.filter.num, tx.filter.den, 2^15, fs); % equivalent two-sided noise bandwidth over larger number of samples (2^15) given a sampling frequency fs 

%% Receiver parameters
% Antialiasing filter
rx.filter = design_filter(rx.filter.type, rx.filterN, 1/(ofdm.Ms*sim.Mct));       

%% define clipping ratios based on optimization results
% this script creates the variables sim.rcliptx, sim.rcliprx, RCLIPTX and
% RCLIPRX based on optimized values and on whether or not quantization is
% on
% [tx.rclip, rx.rclip] = select_clipping_ratio(sim.type, ofdm.CS, tx.modulator.fc);
tx.rclip = sim.rclip;
rx.rclip = sim.rclip;

for k = 1:length(tx.modulator.Fc)
    %% Transmitter parameters 
    tx.modulator.fc = tx.modulator.Fc(k); % modulator cut off frequency
    tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
    tx.modulator.h = @(t) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0));
    tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds                
    
    [ber, Ptx, Ptx_est] = dc_ofdm(ofdm, tx, fiber, rx, sim);
    
    PtxdBm(k) = 10*log10(Ptx/1e-3);
    PtxedBm(k) = 10*log10(Ptx_est/1e-3);
    bercount(k) = log10(ber.count);
    berest(k) = log10(ber.est);
    
    %% Required power for OOK AWGN
    PookdBm(k) = 10*log10(1/rx.R*qfuncinv(sim.BERtarget)*sqrt(rx.Sth*ofdm.Rb)/1e-3) + fiber.att(tx.lamb)*fiber.L/1e3;

    % Same as solving for square constellations
    % qfuncinv((sim.BERtarget*log2(ofdm.CS)/4)/(1-1/sqrt(ofdm.CS)))/sqrt(3*pi*log2(ofdm.CS)/((ofdm.CS-1)*4*rx.Sth*ofdm.Rb))
end

pp.power_pen_ook_m = PtxdBm - PookdBm;
pp.power_pen_ook_e = PtxedBm - PookdBm;

figure, hold on, grid on
plot(tx.modulator.Fc/1e9, bercount, '-sk')
plot(tx.modulator.Fc/1e9, berest, '-or')
plot(tx.modulator.Fc/1e9, log10(sim.BERtarget)*ones(size(tx.modulator.Fc)), 'k')
xlabel('Modulator Cutoff Frequency (GHz)', 'FontSize', 12)
ylabel('log_{10}(BER)', 'FontSize', 12)
legend('BER counted', 'BER estimaded', 'Target BER', 'Location', 'NorthWest')
axis([tx.modulator.Fc([1 end])/1e9 -8 0])
