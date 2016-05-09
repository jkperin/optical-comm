function ber = ber_coherent_analog_dpsk(Tx, Fiber, Rx, sim)
%% Calculate BER of coherent system. Phase estimation is done via electric PLL based on Costas loop
dataTX = randi([0 sim.M-1], [2, sim.Nsymb]); % symbol stream for each polarization

[Vin, symbolsTX] = QAM_SC_Tx(dataTX, Tx, sim); % generates QAM signal

% Modulation format
M = sim.M; % Modulation order
demodulate = @(X) dpskdemod(1/sqrt(2)*exp(-1j*pi/4)*X, M, 0, 'gray').';
% 1/sqrt(2)*exp(-1j*pi/4) is to generate same constellation as if using 4-QAM
if strcmpi(sim.ModFormat, 'QAM')
    error('ber_coherent_analog_dqpsk only supports DPSK');
end

Analog = Rx.Analog;

% BER of ideal system
[ber.theory, ber.theory_noise_enhancement] = ber_coherent_awgn(Tx, Fiber, Rx, sim);

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
          
    %% Pol demux
    %% !!! Must be implemented 

    %% Detection
    % Sampling
    Y = Ys(:, floor(sim.Mct/2):sim.Mct:end);
                   
    % Detection
    % Finds point where DPSK transmission started so that differential
    % decoder starts correctly
    Nstart = find(validInd == sim.Nsetup+1);
    dataRX = [demodulate(Y(1, Nstart:end)); demodulate(Y(2, Nstart:end))];
    validInd(validInd < sim.Nsetup+1) = []; 
    
    dataRX(:, [1:sim.Ndiscard end-sim.Ndiscard+1:end]) = []; % discard first and last sim.Ndiscard symbols
    validInd([1:sim.Ndiscard end-sim.Ndiscard+1:end]) = []; 
    
    % BER calculation
%     figure, stem(dataTX(1, validInd) ~= dataRX(1, :))
%     drawnow
    [~, berX(k)] = biterr(dataTX(1, validInd), dataRX(1, :))
    [~, berY(k)] = biterr(dataTX(2, validInd), dataRX(2, :))
    ber.count(k) = 0.5*(berX(k) + berY(k));

   % Stop simulation when counted BER reaches 0
   if sim.stopWhenBERreaches0 && ber.count(k) == 0
       break
   end
end

plots

%% Plot BER curves
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
    
