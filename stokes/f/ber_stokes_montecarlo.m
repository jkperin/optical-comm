function ber = ber_stokes_montecarlo(Tx, Fiber, Rx, sim)
%% Calculate BER of coherent system with DSP-based receiver for values of launched power specified Tx.PlaunchdB
ModFormat = sim.ModFormat;

%% Generate transmitted symbols
% Note 1: this implement assumes that no DAC is being used so ModFormat.signal
% generates rectangular pulses sampled at sim.fs
% Note 2: ModFormat.signal does not remove group delay due to pulse shaping
% so the center of the pulse is at (Mct+1)/2, if pulse shape is 'rect'
dataTX = randi([0 ModFormat(1).M-1], [2, sim.Nsymb]); % symbol stream for each polarization
Nzero = 10; % zero Nzero first and last symbols to make sequence periodic
dataTX(:, 1:Nzero) = 0;
dataTX(:, end-Nzero+1:end) = 0;

ModFormat = ModFormat.unbias;
[Vin(1, :), symbolsTX(1, :)] = ModFormat.signal(dataTX(1, :)); % X & Y pol
[Vin(2, :), symbolsTX(2, :)] = ModFormat.signal(dataTX(2, :)); % X & Y pol

% Filter drive waveforms for modulators txfilt.
% group delay of Tx.filt.H has already been removed
Vout(1, :) = dac(Vin(1, :), Tx.DAC, sim);
Vout(2, :) = dac(Vin(2, :), Tx.DAC, sim);

%% Swipe launched power
% Transmitted power swipe
ber.count = zeros(size(Tx.PlaunchdBm));
for k = 1:length(Tx.PlaunchdBm)
    fprintf('-- Launch power: %.2f dBm\n', Tx.PlaunchdBm(k));
    
    %% ========= Modulator ==========
    Tx.Laser.PdBm = Tx.PlaunchdBm(k);
    Ecw = Tx.Laser.cw(sim); % Generates electric field with intensity and phase noise   
    Vout = Tx.Mod.Vswing*(Vout/2) + Tx.Mod.Vbias;
    Ein(1, :) = mzm(Ecw/sqrt(2), Vout(1, :), Tx.Mod); % modulate optical signal using eletro-optical modulator (EOM)
    Ein(2, :) = mzm(Ecw/sqrt(2), Vout(2, :), Tx.Mod); % modulate optical signal using eletro-optical modulator (EOM)
    
    % Ensure that transmitted power is a desired level
    Ein = Ein*sqrt(dBm2Watt(Tx.PlaunchdBm(k))/sum(mean(abs(Ein).^2, 2)));

    %% ========= Propagation ========== 
    Erec = Fiber.linear_propagation(Ein, sim.f, Tx.Laser.lambda);
    
    if isfield(sim, 'preAmp') && sim.preAmp % only included if sim.preAmp is true
        disp('- IMPORTANT: Simulation including optical amplifier!')
        [Erec, OSNRdBtheory] = Rx.OptAmp.amp(Erec, sim.fs);

        % Measure OSNR
        Osa = OSA(0.1); % optical spectrum analyser with resolution 0.1nm
        OSNRdBmeasured = Osa.estimate_osnr(Erec, Tx.Laser.wavelength, sim.f, sim.shouldPlot('OSNR'));

        fprintf('OSNR = %.2f dB (theory)\nOSNR = %.2f dB (measured)\n', OSNRdBtheory, OSNRdBmeasured)
               
        % check
%         OSNRdBtheorycheck = 10*log10(dBm2Watt(Tx.PlaunchdBm(k))/(2*Rx.OptAmp.nsp*Rx.OptAmp.h*Rx.OptAmp.c/Tx.Laser.lambda*12.5e9))
%         SNRdBtheorycheck = 10*log10(dBm2Watt(Tx.PlaunchdBm(k))/(Rx.OptAmp.nsp*Rx.OptAmp.h*Rx.OptAmp.c/Tx.Laser.lambda*sim.Rs))
    end
    
    %% ========= Receiver =============
    Yrx = dual_pol_stokes_receiver(Erec, Rx, sim);
    
    %% Sampling    
    % Digital to analog conversion: antialiasing filter, sampling, and
    % quantization
    Y(1, :) = adc(Yrx(1, :), Rx.ADC, sim);
    Y(2, :) = adc(Yrx(2, :), Rx.ADC, sim);
    Y(3, :) = adc(Yrx(3, :), Rx.ADC, sim);
    Y(4, :) = adc(Yrx(4, :), Rx.ADC, sim);
    
    %% Automatic gain control
    Y = Y*sqrt(mean(abs(symbolsTX(:)).^2)/mean(abs(Y(:)).^2));
       
    %% Equalization
    Rx.eq.trainSeq = symbolsTX;
    X = equalize_mimo4x2(Y, Rx.eq, ModFormat);
    
    %% Detection
    dataRX(1, :) = ModFormat.demod(X(1, sim.Ndiscard+Rx.eq.Ndiscard+1:end-sim.Ndiscard));
    dataRX(2, :) = ModFormat.demod(X(2, sim.Ndiscard+Rx.eq.Ndiscard+1:end-sim.Ndiscard));

    %% BER calculation
    validRange = sim.Ndiscard+Rx.eq.Ndiscard+1:sim.Nsymb-sim.Ndiscard;
    [~, berX(k)] = biterr(dataTX(1, validRange), dataRX(1, :))
    [~, berY(k)] = biterr(dataTX(2, validRange), dataRX(2, :))
    ber.count(k) = 0.5*(berX(k) + berY(k));
   
   % Stop simulation when counted BER reaches 0
   if isfield(sim, 'stopWhenBERreaches0') && sim.stopWhenBERreaches0 && ber.count(k) == 0
       break
   end   
end

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
    
