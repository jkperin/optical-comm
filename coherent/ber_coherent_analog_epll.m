function ber = ber_coherent_analog_epll(Tx, Fiber, Rx, sim)
%% Calculate BER of coherent system. Phase estimation is done via electric PLL based on XOR operations
dataTX = randi([0 sim.M-1], [2, sim.Nsymb]); % symbol stream for each polarization

[Vin, symbolsTX] = QAM_SC_Tx(dataTX, Tx, sim); % generates QAM signal

% Modulation format
M = sim.M; % Modulation order
demodulate = @(X) qamdemod(X, M, 0, 'gray');
if strcmpi(sim.ModFormat, 'DPSK') || M ~= 4
    error('ber_coherent_analog_epll_costas only supports QPSK');
end

% BER of ideal system
[ber.theory, ber.theory_noise_enhancement] = ber_coherent_awgn(Tx, Fiber, Rx, sim);

Analog = Rx.Analog;

% Transmitted power swipe
ber.count = zeros(size(Tx.PlaunchdBm));
for k = 1:length(Tx.PlaunchdBm)
    fprintf('Launch power: %.2f dBm\n', Tx.PlaunchdBm(k));
    validInd = 1:sim.Nsymb;
    
    Tx.Laser.PdBm = Tx.PlaunchdBm(k);
    
    Ein = Tx.Laser.cw(sim); % Generates electric field with intensity and phase noise
       
    if strcmpi(sim.Modulator, 'MZM')
        Ein = mzm(Ein, Vin, Tx.Mod); % modulate optical signal using eletro-optical modulator (EOM)
    else
        Ein = SiPh_optical_mod(Ein, Vin, Tx.Mod);
    end
    
    % Ensure that transmitted power is a desired level
    Ein = Ein*sqrt(dBm2Watt(Tx.PlaunchdBm(k))/sum(mean(abs(Ein).^2, 2)));

    %% ========= Propagation ========== 
    Erec = Fiber.linear_propagation(Ein, sim.f, Tx.Laser.lambda);
    
    %% ========= Receiver =============
    Y = dual_pol_coherent_receiver(Erec, sim.M, Rx, sim);
    
    % Receiver filter
    Hrx = ifftshift(Analog.filt.H(sim.f/sim.fs));
    Yxi = real(ifft(fft(real(Y(1, :))).*Hrx));
    Yxq = real(ifft(fft(imag(Y(1, :))).*Hrx));
    Yyi = real(ifft(fft(real(Y(2, :))).*Hrx));
    Yyq = real(ifft(fft(imag(Y(2, :))).*Hrx));

    Ys = [Yxi + 1j*Yxq; Yyi + 1j*Yyq];
          
    %% Pol demux !!! Must be implemented 
    
    %% Carrier phase recovery via electric PLL (EPLL)
    switch lower(Analog.CPRmethod)
        case 'costas' % EPLL via Costas loop
           [Xs, LoopFilterInput, LoopFilterOutput, Analog] = analog_epll_costas(Ys, Tx.Laser.linewidth + Rx.LO.linewidth, Analog, sim);
        case 'logic' % EPLL via logic (XOR) operations
           [Xs, LoopFilterInput, LoopFilterOutput, Analog] = analog_epll_logic(Ys, Tx.Laser.linewidth + Rx.LO.linewidth, Analog, sim);
        otherwise
            error('ber_coherent_analog_epll/invalid electric PLL type %s\nAnalog.receiver must be either Costas or Logic\n', Analog.receiver)
    end
    
    % Sampling
    Y = Xs(:, sim.Mct/2:sim.Mct:end);  
    % Note: decisions could be obtained directly from Xd. Xs is used here 
    % instead so we can see the constellation 
    
    % Automatic gain control
    Y(1, :) = sqrt(2/mean(abs(Y(1, :)).^2))*Y(1, :);
    Y(2, :) = sqrt(2/mean(abs(Y(2, :)).^2))*Y(2, :);
    
    % Correct for constant phase rotation
    phaseShifts = [0 pi/2 pi 3*pi/2 pi/4+[0 pi/2 pi 3*pi/2]];
    Npoints = sim.Ndiscard:sim.Ndiscard+1e3;
    errors = zeros(2, length(phaseShifts));
    for ind = 1:length(phaseShifts)
        errors(1, ind) = biterr(demodulate(Y(1, Npoints)*exp(1j*phaseShifts(ind))), dataTX(1, Npoints));
        errors(2, ind) = biterr(demodulate(Y(2, Npoints)*exp(1j*phaseShifts(ind))), dataTX(2, Npoints));
    end
     
    [~, ind1] = min(errors(1, :));
    [~, ind2] = min(errors(2, :));
    
    Y(1, :) = exp(1j*phaseShifts(ind1))*Y(1, :);
    Y(2, :) = exp(1j*phaseShifts(ind2))*Y(2, :);
    
    % Detection
    dataRX = [demodulate(Y(1, :)); demodulate(Y(2, :))];
    
    dataRX(:, [1:sim.Ndiscard end-sim.Ndiscard+1:end]) = []; % discard first and last sim.Ndiscard symbols
    validInd([1:sim.Ndiscard end-sim.Ndiscard+1:end]) = []; 
    
    % BER calculation
%     figure, stem(dataTX(1, validInd) ~= dataRX(1, :))
%     drawnow
    [~, berX(k)] = biterr(dataTX(1, validInd), dataRX(1, :))
    [~, berY(k)] = biterr(dataTX(2, validInd), dataRX(2, :))
    ber.count(k) = 0.5*(berX(k) + berY(k));
    
   % Constellation plots
   if sim.Plots.isKey('Constellations') && sim.Plots('Constellations')
       figure(203), clf
       subplot(121)
       plot_constellation(Y(1, sim.Ndiscard+1:end-sim.Ndiscard),...
           dataTX(1, sim.Ndiscard+1:end-sim.Ndiscard), sim.M);
       axis square
       title('Pol X')   
      subplot(122)
       plot_constellation(Y(2, sim.Ndiscard+1:end-sim.Ndiscard),...
           dataTX(2, sim.Ndiscard+1:end-sim.Ndiscard), sim.M);
       axis square
       title('Pol Y')   
       drawnow
   end 
   
   % Stop simulation when counted BER reaches 0
   if sim.stopWhenBERreaches0 && ber.count(k) == 0
       break
   end
end

plots

%% Plot BER curve
if isfield(sim, 'Plots') && sim.Plots('BER') && length(ber.count) > 1
    figure(100), box on, hold on
    [~, link_attdB] = Fiber.link_attenuation(Tx.Laser.lambda);
    Prx = Tx.PlaunchdBm - link_attdB;
    plot(Prx, log10(ber.theory), '-')
    plot(Prx, log10(ber.count), '-o')
    legend('Theory', 'Counted')
    xlabel('Received Power (dBm)')
    ylabel('log_{10}(BER)')
    axis([Prx(1) Prx(end) -8 0])
end
    
