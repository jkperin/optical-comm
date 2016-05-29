function ber = ber_coherent_dsp(Tx, Fiber, Rx, sim)
%% Calculate BER of coherent system
dataTX = randi([0 sim.M-1], [2, sim.Nsymb]); % symbol stream for each polarization

[Vin, symbolsTX] = QAM_SC_Tx(dataTX, Tx, sim); % generates QAM signal

% Modulation format
M = sim.M; % Modulation order
if strcmpi(sim.ModFormat, 'QAM')
    demodulate = @(X) qamdemod(X, M, 0, 'gray');
elseif strcmpi(sim.ModFormat, 'DPSK')
    demodulate = @(X) dpskdemod(1/sqrt(2)*exp(-1j*pi/4)*X, M, 0, 'gray').';
    % 1/sqrt(2)*exp(-1j*pi/4) is to generate same constellation as if using 4-QAM
else
    error('ber_coherent/invalid modulation format');
end

% BER of ideal system
[ber.theory, ~, SNRdBtheory] = ber_coherent_awgn(Tx, Fiber, Rx, sim);

% Transmitted power swipe
ber.count = zeros(size(Tx.PlaunchdBm));
for k = 1:length(Tx.PlaunchdBm)
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
    
    % Digital to analog conversion: antialiasing filter, sampling, and
    % quantization
    [Yxi, varQ(1)] = adc(real(Y(1, :)), Rx.ADC, sim);
    [Yxq, varQ(2)] = adc(imag(Y(1, :)), Rx.ADC, sim);
    [Yyi, varQ(3)] = adc(real(Y(2, :)), Rx.ADC, sim);
    [Yyq, varQ(4)] = adc(imag(Y(2, :)), Rx.ADC, sim);

    Ys = [Yxi + 1j*Yxq; Yyi + 1j*Yyq];
          
    %% Equalization
    switch upper(Rx.AdEq.type) 
        case 'CMA'
            [Yeq, W, MSE] = fse_cma(Ys, Rx.AdEq, sim); % equalization
        case 'LMS'
            [Yeq, W, MSE] = lms_td_fse(Ys, Rx.AdEq, symbolsTX, sim.M);
        otherwise
            error('ber_coherent/Unknown Rx.AdEq.type option')
    end
    validInd(1:Rx.AdEq.Ntrain) = []; % discard symbols used in training equalizer
    Yeq(:, 1:Rx.AdEq.Ntrain) = [];    

    %% Carrier phase recovery algorithm
    CPRsetup = 0;
    if strcmpi(sim.ModFormat, 'QAM') % not done if DQPSK
        Rx.CPR.Ytrain = symbolsTX(:, validInd);
        if strcmpi(Rx.CPR.type, 'feedforward')
            Rx.CPR.SNRdB = SNRdBtheory(k);
            Rx.CPR.varPN = Tx.Laser.varPN(sim.Rs) + Rx.LO.varPN(sim.Rs);
            
            % Frequency estimation
            Ycpr = Yeq;
            if Rx.LO.freqOffset ~= 0
                [Ycpr, Foff_est] = qpsk_freq_rec(Ycpr, Rx.FreqRec, sim);
                Ycpr(:, 1:Rx.FreqRec.Ntrain) = []; % discard symbols used in training phase tracking algorithm     
                validInd(1:Rx.FreqRec.Ntrain) = [];         
                CPRsetup = CPRsetup + Rx.FreqRec.Ntrain;
            end
            
            % Feedforward carrier recovery
            Ycpr = feedforward_cpr(Ycpr, Rx.CPR, sim);      
            Ycpr(:, 1:Rx.CPR.Ntrain) = [];
            Y = Ycpr;
            validInd(1:Rx.CPR.Ntrain) = []; % discard symbols used in training carrier phase recovery algorithm
            CPRsetup = CPRsetup + Rx.CPR.Ntrain;
            
            % Phase tracking (constant phase offset compensation)
            [Y, phiPT] = phase_estimation(Y, Rx.PT, symbolsTX(:, validInd), sim);
            Y(:, 1:Rx.PT.Ntrain) = [];
            validInd(1:Rx.PT.Ntrain) = []; % discard symbols used in training phase tracking algorithm            
        elseif strcmpi(Rx.CPR.type, 'DPLL')
            % DPLL
            Ycpr = dpll(Yeq, Rx.CPR, sim);
            Ycpr(:, 1:Rx.CPR.Ntrain) = [];
            Y = Ycpr;
            validInd(1:Rx.CPR.Ntrain) = []; % discard symbols used in training carrier phase recovery algorithm
            CPRsetup = CPRsetup + Rx.CPR.Ntrain;
            
            % Phase tracking (constant phase offset compensation)
            [Y, phiPT] = phase_estimation(Y, Rx.PT, symbolsTX(:, validInd), sim);
            Y(:, 1:Rx.PT.Ntrain) = [];
            validInd(1:Rx.PT.Ntrain) = []; % discard symbols used in training phase tracking algorithm          
        elseif strcmpi(Rx.CPR.type, 'None')
            Rx.CPR.Ntrain = 0;
            Ycpr = Yeq;
            Y = Yeq;
        else
            error('ber_coherent/invalid carrier phase recovery type')
        end        
    elseif strcmpi(sim.ModFormat, 'DPSK')
        Y = Yeq;
        if Rx.LO.freqOffset ~= 0
            [Y, Foff_est] = qpsk_freq_rec(Y, Rx.FreqRec, sim);
            Y(:, 1:Rx.FreqRec.Ntrain) = []; % discard symbols used in training phase tracking algorithm     
            validInd(1:Rx.FreqRec.Ntrain) = []; 
            CPRsetup = CPRsetup + Rx.FreqRec.Ntrain;
        end
    end
    
    % Automatic gain control
    Y(1, :) = sqrt(2/mean(abs(Y(1, :)).^2))*Y(1, :);
    Y(2, :) = sqrt(2/mean(abs(Y(2, :)).^2))*Y(2, :); 
    
    % Detection
    if strcmpi(sim.ModFormat, 'QAM')
        dataRX = [demodulate(Y(1, :)); demodulate(Y(2, :))];
    elseif strcmpi(sim.ModFormat, 'DPSK')
        % Finds point where DPSK transmission started so that differential
        % decoder starts correctly
        Nstart = find(validInd == sim.Nsetup+1);
        dataRX = [demodulate(Y(1, Nstart:end)); demodulate(Y(2, Nstart:end))];
        validInd(validInd < sim.Nsetup+1) = []; 
    end
    
    dataRX(:, [1:sim.Ndiscard end-sim.Ndiscard+1:end]) = []; % discard first and last sim.Ndiscard symbols
    validInd([1:sim.Ndiscard end-sim.Ndiscard+1:end]) = []; 
    
    % BER calculation
%     figure, stem(dataTX(1, validInd) ~= dataRX(1, :))
%     drawnow
    [~, berX(k)] = biterr(dataTX(1, validInd), dataRX(1, :))
    [~, berY(k)] = biterr(dataTX(2, validInd), dataRX(2, :))
    ber.count(k) = 0.5*(berX(k) + berY(k));

   % Constellation plots
   if strcmpi(sim.ModFormat, 'QAM') && sim.Plots.isKey('Constellations') && sim.Plots('Constellations')
       figure(203), clf
       subplot(121)
       plot_constellation(Yeq(1, sim.Ndiscard+1:end-sim.Ndiscard),...
           dataTX(1, Rx.AdEq.Ntrain+sim.Ndiscard+1:end-sim.Ndiscard), sim.M);
       axis square
       title('Pol X: After equalization')
       subplot(122)
       plot_constellation(Ycpr(1, sim.Ndiscard+1:end-sim.Ndiscard),...
           dataTX(1, Rx.AdEq.Ntrain+CPRsetup+sim.Ndiscard+1:end-sim.Ndiscard), sim.M);
       axis square
       title('Pol X: After CPR')   
       drawnow
   end 
   
   % Stop simulation when counted BER reaches 0
   if isfield(sim, 'stopWhenBERreaches0') && sim.stopWhenBERreaches0 && ber.count(k) == 0
       break
   end   
end

plots

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
    
