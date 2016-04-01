function [ber, SNRdB_NE] = ber_coherent(Tx, Fiber, Rx, sim)
%% Simulate transmission of cohrent system
% sequenceX = debruijn_sequence(sim.M, 5);
% sequenceY = sequenceX(end:-1:1);
% sequence = repmat([sequenceX; sequenceY], 1, 4);
dataTX = randi([0 sim.M-1], [2, sim.Nsymb]); % symbol stream for each polarization

[Vin, symbolsTX] = QAM_SC_Tx(dataTX, Tx, sim); % generates QAM signal

% Modulation format
M = sim.M; % Modulation order
if strcmpi(sim.ModFormat, 'QAM');
    demodulate = @(X) qamdemod(X, M, 0, 'gray');
elseif strcmpi(sim.ModFormat, 'DPSK');
    demodulate = @(X) dpskdemod(1/sqrt(2)*exp(-1j*pi/4)*X, M, 0, 'gray').';
    % 1/sqrt(2)*exp(-1j*pi/4) is to generate same constellation as if using 4-QAM
else
    error('ber_coherent/invalid modulation format');
end

% BER in AWGN channel
berAWGN = @(SNRdB) berawgn(SNRdB - 10*log10(log2(M)), lower(sim.ModFormat), M);
Df = (Rx.ADC.filt.noisebw(sim.fs)/2)/(Rx.ADC.fs/2);

% Transmitted power swipe
for k = 1:length(Tx.PlaunchdBm)
    validInd = 1:sim.Nsymb;
    
    Tx.Laser.PdBm = Tx.PlaunchdBm(k);
    
    Ein = Tx.Laser.cw(sim); % Generates electric field with intensity and phase noise

    if strcmpi(sim.Modulator, 'EOM')
        Ein = eom(Ein, Vin, Tx.Mod); % modulate optical signal using eletro-optical modulator (EOM)
    else
        Ein = SiPh_optical_mod(Ein, Vin, Tx.Mod);
    end

    %% ========= Propagation ========== 
    Erec = Fiber.linear_propagation(Ein, sim.f, Tx.Laser.lambda);

    %% ========= Receiver =============
    [Y, ~] = PDM_QAM_Rx(Erec, sim.M, Rx, sim);
    
    % Digital to analog conversion: filters, sample, and quantize
    [Yxi, varQ(1)] = adc(real(Y(1, :)), Rx.ADC, sim);
    [Yxq, varQ(2)] = adc(imag(Y(1, :)), Rx.ADC, sim);
    [Yyi, varQ(3)] = adc(real(Y(2, :)), Rx.ADC, sim);
    [Yyq, varQ(4)] = adc(imag(Y(2, :)), Rx.ADC, sim);

    Ys = [Yxi + 1j*Yxq; Yyi + 1j*Yyq];
    
    Prx = dBm2Watt(Tx.Laser.PdBm)/Fiber.link_attenuation(Tx.Laser.lambda);
    Plo = dBm2Watt(Rx.LO.PdBm);
    
    %% Equalization
    switch upper(Rx.AdEq.type) 
        case 'CMA'
            [Y, W, MSE] = fse_cma(Ys, Rx.AdEq, sim); % equalization
        case 'LMS'
            [Y, W, MSE] = lms_td_fse(Ys, Rx.AdEq, symbolsTX, sim.M);
        otherwise
            error('ber_coherent/Unknown Rx.AdEq.type option')
    end
    validInd(1:Rx.AdEq.Ntrain) = []; % discard symbols used in training equalizer
    
    % Estimate SNR including noise enhacement penalty
%     noiseBW = Rx.ADC.fs*calc_noiseBW(W{1}, W{2}); % mean noise bandwidth of equalizers
    noiseBW = sim.Rs/2; % noise bandwidth if receiver filter were matched filter

    Ppd = abs(sqrt(Plo/(4*sim.Npol)) + sqrt(Prx/(4*sim.Npol))).^2; % incident power in each photodiode
    Psig = 1/sqrt(2)*Plo*Prx/(sim.Npol*sim.Npol); % Signal power per polarization
    varShot = 4*Rx.PD.varShot(Ppd, noiseBW*Df);
    varThermal = 2*Rx.N0*noiseBW*Df;
    SNRdB_NE(k) = 10*log10(Psig/(varShot + varThermal));
           
    %% Carrier phase recovery algorithm
    if isfield(Rx, 'CPR') && strcmpi(sim.ModFormat, 'QAM') % not done if DQPSK
        Rx.CPR.Ytrain = symbolsTX(:, validInd);
        if strcmpi(Rx.CPR.type, 'feedforward')
            Rx.CPR.SNRdB = SNRdB_NE(k);
            Rx.CPR.varPN = Tx.Laser.varPN(sim.Rs) + Rx.LO.varPN(sim.Rs);
            Yp = feedforward_cpr(Y(:, validInd), Rx.CPR, sim);      
        elseif strcmpi(Rx.CPR.type, 'DPLL')
            Yp = dpll(Y(:, validInd), Rx.CPR, sim);
        else
            error('ber_coherent/invalid carrier phase recovery type')
        end
        Y = [Y(:, 1:Rx.AdEq.Ntrain) Yp];
        validInd(1:Rx.CPR.Ntrain) = []; % discard symbols used in training carrier phase recovery algorithm
        
        % Phase tracking
        [Yp, phiPT] = phase_estimation(Y(:, validInd), Rx.PT, symbolsTX(:, validInd), sim);
        Y = [Y(:, 1:Rx.AdEq.Ntrain+Rx.CPR.Ntrain) Yp];
        validInd(1:Rx.PT.Ntrain) = []; % discard symbols used in training phase tracking algorithm
    elseif strcmpi(sim.ModFormat, 'DPSK')
        [Yp, Foff_est] = qpsk_freq_rec(Y(:, validInd), Rx.FreqRec);
        mean(Foff_est)*sim.Rs
        Y = [Y(:, 1:Rx.AdEq.Ntrain) Yp];
        validInd(1:Rx.FreqRec.Ntrain) = []; % discard symbols used in training phase tracking algorithm     
    end
    
    validInd([1:sim.Ndiscard end-sim.Ndiscard+1:end]) = []; % discard first and last sim.Ndiscard symbols
    dataRX = [demodulate(Y(1, :)); demodulate(Y(2, :))];
    [~, berX(k)] = biterr(dataTX(1, validInd), dataRX(1, validInd))
    [~, berY(k)] = biterr(dataTX(2, validInd), dataRX(2, validInd))
    ber(k) = 0.5*(berX(k) + berY(k));
    berAWGN(SNRdB_NE(k));
end

plots

if isfield(sim, 'Plots') && sim.Plots('BER') && length(ber) > 1
    figure(1), box on, hold on
    [~, link_attdB] = Fiber.link_attenuation(Tx.Laser.lambda);
    Prx = Tx.PlaunchdBm - link_attdB;
    plot(Prx, log10(berawgn(SNRdB_NE - 10*log10(log2(sim.M)), lower(sim.ModFormat), sim.M)), '-')
    plot(Prx, log10(ber), '-o')
    legend('Theory including noise enhancement', 'Counted')
    xlabel('Received Power (dBm)')
    ylabel('log_{10}(BER)')
    axis([Prx(1) Prx(end) -8 0])
end
    
