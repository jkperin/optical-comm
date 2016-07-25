function ber = ber_coherent_analog_dpsk(Tx, Fiber, Rx, sim)
%% Calculate BER of DP-DQPSK coherent system with analog-based receiver.

Dpsk = sim.ModFormat;

% Ensures that modulation format is DPSK
assert(strcmpi(class(Dpsk), 'DPSK'), 'ber_coherent_analog_dpsk: Modulation format must be DPSK')

%% Generate transmitted symbols
% Note 1: this implement assumes that no DAC is being used so ModFormat.signal
% generates rectangular pulses sampled at sim.fs
% Note 2: ModFormat.signal does not remove group delay due to pulse shaping
% so the center of the pulse is at (Mct+1)/2, if pulse shape is 'rect'
dataTX = randi([0 Dpsk.M-1], [2, sim.Nsymb]); % symbol stream for each polarization
Nzero = 10; % zero Nzero first and last symbols to make sequence periodic
dataTX(:, 1:Nzero) = 0;
dataTX(:, end-Nzero+1:end) = 0;

[Vin, symbolsTX] = Dpsk.signal(dataTX); % X & Y pol

% Filter drive waveforms for modulators txfilt.
% group delay of Tx.filt.H has already been removed
Htx = ifftshift(Tx.filt.H(sim.f/sim.fs).*exp(1j*2*pi*sim.f/sim.fs*(Dpsk.pulse_shape_grpdelay))); % transmitter filter and remove group delay due to pulse shaping in ModFormat
Vout(1, :) = real(ifft(fft(real(Vin(1, :))).*Htx)) + 1j*real(ifft(fft(imag(Vin(1, :))).*Htx)); 
Vout(2, :)= real(ifft(fft(real(Vin(2, :))).*Htx)) + 1j*real(ifft(fft(imag(Vin(2, :))).*Htx));

% BER of ideal system with receiver filtering equal to matched filter
[ber.theory, ber.theory_noise_enhancement] = ber_coherent_awgn(Tx, Fiber, Rx, sim);

%% Swipe launched power
ber.count = zeros(size(Tx.PlaunchdBm));
for k = 1:length(Tx.PlaunchdBm)
    fprintf('-- Launch power: %.2f dBm\n', Tx.PlaunchdBm(k));
    
    %% ========= Modulator ==========
    Tx.Laser.PdBm = Tx.PlaunchdBm(k);
    Ein = Tx.Laser.cw(sim); % Generates electric field with intensity and phase noise   
    Ein = mzm(Ein, Vout, Tx.Mod); % modulate optical signal using eletro-optical modulator (EOM)
    
    % Ensure that transmitted power is a desired level
    Ein = Ein*sqrt(dBm2Watt(Tx.PlaunchdBm(k))/sum(mean(abs(Ein).^2, 2)));

    %% ========= Propagation ========== 
    Fiber.PMD = false; 
    % Note: since polarization demultiplexing is not done here, the fiber
    % must maitain polarization states.
    Erec = Fiber.linear_propagation(Ein, sim.f, Tx.Laser.lambda);
    
    %% ========= Receiver =============
    Yrx = dual_pol_coherent_receiver(Erec, Dpsk.M, Rx, sim);
    
    %% Receiver filter
    Hrx = ifftshift(Rx.Analog.filt.H(sim.f/sim.fs));
    Ys = [real(ifft(fft(real(Yrx(1, :))).*Hrx)) + 1j*real(ifft(fft(imag(Yrx(1, :))).*Hrx));...
          real(ifft(fft(real(Yrx(2, :))).*Hrx)) + 1j*real(ifft(fft(imag(Yrx(2, :))).*Hrx))];
    

    %% Time recovery and sampling
    % Note: clock obtained from I is used in Q in order to prevent
    % differences in sampling between I and Q
    [X(1, :), ~, Nsetup] = analog_time_recovery(Ys(1, :), Rx.TimeRec, sim, sim.shouldPlot('Time recovery'));
    X(2, :) = analog_time_recovery(Ys(2, :), Rx.TimeRec, sim);
       
    % Automatic gain control
    X = [sqrt(2/mean(abs(X(1, :)).^2))*X(1, :);
         sqrt(2/mean(abs(X(2, :)).^2))*X(2, :)];

    %% Align received sequence and correct for phase rotation
    [c(1, :), ind] = xcorr(symbolsTX(1, :), X(1, :), 20, 'coeff');
    c(2, :) = xcorr(symbolsTX(2, :), X(2, :), 20, 'coeff');
       
    % maximum correlation position
    [~, p] = max(abs(c), [], 2);
    theta = [angle(c(1, p(1))), angle(c(2, p(2)))];
    
    % Circularly shift symbol sequence
    X = [circshift(X(1, :), [0 ind(p(1))]);...
        circshift(X(2, :), [0 ind(p(2))])];
       
    % Rotate constellations
    X = [X(1, :).*exp(+1j*theta(1)); X(2, :).*exp(+1j*theta(2))];
    
    %% Detection
    % Note: since DPSK has differential enconding, detection for BER
    % measurement is done after several symbols so that measurements are
    % not currepted by adaptitation, fft windowing, etc
    
    % Valid range for BER measurement
    validInd = sim.Ndiscard+Nsetup(1)+1:sim.Nsymb-sim.Ndiscard-Nsetup(2);
    
    dataRX = Dpsk.demod(X(:, validInd), symbolsTX(:, validInd(1)-1));
        
    % BER calculation
    [~, berX(k)] = biterr(dataTX(1, validInd), dataRX(1, :))
    [~, berY(k)] = biterr(dataTX(2, validInd), dataRX(2, :))
    ber.count(k) = 0.5*(berX(k) + berY(k));

   if sim.shouldPlot('Symbol errors')
        figure(204), clf
        subplot(121)
        stem(dataTX(1, validInd) ~= dataRX(1, validInd))
        title('Pol X')       
        subplot(122)
        stem(dataTX(2, validInd) ~= dataRX(2, validInd))
        title('Pol Y')   
        drawnow
   end
    
   % Stop simulation when counted BER reaches 0
   if sim.stopWhenBERreaches0 && ber.count(k) == 0
       break
   end
end

plots

%% Plot BER curves
if sim.shouldPlot('BER') && length(ber.count) > 1
    figure(100), box on, hold on
    [~, link_attdB] = Fiber.link_attenuation(Tx.Laser.lambda);
    Prx = Tx.PlaunchdBm - link_attdB;
    hline(1) = plot(Prx, log10(ber.theory), '-');
    hline(2) = plot(Prx, log10(ber.count), '-o');
    legend(hline, {'Theory', 'Counted'})
    xlabel('Received Power (dBm)')
    ylabel('log_{10}(BER)')
    axis([Prx(1) Prx(end) -8 0])
end
    
