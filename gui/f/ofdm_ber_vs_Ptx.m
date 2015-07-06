function ber = ofdm_ber_vs_Ptx(ofdm, tx, fiber, rx, sim)

%% Transmitter parameters 
tx.modulator.fc = tx.modulator.Fc(1); % modulator cut off frequency
tx.modulator.H = @(f) 1./(1 + 2*1j*f/tx.modulator.fc - (f/tx.modulator.fc).^2);  % laser freq. resp. (unitless) f is frequency vector (Hz)
tx.modulator.h = @(t) (2*pi*tx.modulator.fc)^2*t(t >= 0).*exp(-2*pi*tx.modulator.fc*t(t >= 0));
tx.modulator.grpdelay = 2/(2*pi*tx.modulator.fc);  % group delay of second-order filter in seconds

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

Ptx = 1e-3*10.^(tx.PtxdBm/10);

% Initiliaze variables
bercount = zeros(size(Ptx));
berest = zeros(size(Ptx));

for k = 1:length(Ptx)
    tx.Ptx = Ptx(k);

    % generate, propagate, and detect
    ber = dc_ofdm(ofdm, tx, fiber, rx, sim);

    bercount(k) = ber.count;
    berest(k) = ber.est;
end

ber.count = bercount;
ber.est = berest;