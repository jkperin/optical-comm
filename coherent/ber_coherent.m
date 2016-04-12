function [ber, SNRdB_theory] = ber_coherent(Tx, Fiber, Rx, sim)
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

% BER in AWGN channel
berAWGN = @(SNRdB) berawgn(SNRdB - 10*log10(log2(M)), lower(sim.ModFormat), M);

% Calculate noise bandwidth of an ideal receiver consisting of a matched 
% filter and symbol-rate linear equalization. This is only used to
% estiamate theoretical BER
eq.type = 'fixed td-sr-le';
eq.Ntaps = 31;
mpam = PAM(2, sim.Rs, 'equally-spaced', @(n) double(n >= 0 & n < sim.Mct));
Hch = (Tx.filt.H(sim.f/sim.fs).*Tx.Mod.Hel).';
sim.f = sim.f.';
[~, eq] = equalize(eq, [], Hch, mpam, Rx, sim); 
% this generates zero-forcing linear equalizer, which should be a good
% approximation of MMSE linear equalizer.
noiseBW = trapz(sim.f, abs(eq.Hrx.*eq.Hff(sim.f/sim.Rs)).^2)/2;
sim.f = sim.f.';

% Transmitted power swipe
ber.count = zeros(size(Tx.PlaunchdBm));
ber.theory = zeros(size(Tx.PlaunchdBm));
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
    [Y, ~] = PDM_QAM_Rx(Erec, sim.M, Rx, sim);
    
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
    
    % Estimate SNR including noise enhacement penalty of an ideal symbol-rate linear equalizer   
    Prx = dBm2Watt(Tx.Laser.PdBm)/Fiber.link_attenuation(Tx.Laser.lambda); % received power
    Plo = dBm2Watt(Rx.LO.PdBm);                                            % local oscillator power
    Ppd = abs(sqrt(Plo/(4*sim.Npol)) + sqrt(Prx/(4*sim.Npol))).^2;         % incident power in each photodiode
    Psig = 1/sqrt(2)*Plo*Prx/(sim.Npol*sim.Npol);                          % Signal power per real dimension
    varShot = 2*Rx.PD.varShot(Ppd, noiseBW);                               % Shot noise variance per real dimension
    varThermal = Rx.N0*noiseBW;                                            % Thermal noise variance per real dimension
    SNRdB_theory(k) = 10*log10(Psig/(varShot + varThermal));

    %% Carrier phase recovery algorithm
    if strcmpi(sim.ModFormat, 'QAM') % not done if DQPSK
        Rx.CPR.Ytrain = symbolsTX(:, validInd);
        if strcmpi(Rx.CPR.type, 'feedforward')
            Rx.CPR.SNRdB = SNRdB_theory(k);
            Rx.CPR.varPN = Tx.Laser.varPN(sim.Rs) + Rx.LO.varPN(sim.Rs);
            Ycpr = feedforward_cpr(Yeq, Rx.CPR, sim);      
            
            Ycpr(:, 1:Rx.CPR.Ntrain) = [];
            Y = Ycpr;
            validInd(1:Rx.CPR.Ntrain) = []; % discard symbols used in training carrier phase recovery algorithm
            
            % Phase tracking
            [Y, phiPT] = phase_estimation(Y, Rx.PT, symbolsTX(:, validInd), sim);
            Y(:, 1:Rx.PT.Ntrain) = [];
            validInd(1:Rx.PT.Ntrain) = []; % discard symbols used in training phase tracking algorithm
            
        elseif strcmpi(Rx.CPR.type, 'DPLL')
            Ycpr = dpll(Yeq, Rx.CPR, sim);
            
            Ycpr(:, 1:Rx.CPR.Ntrain) = [];
            Y = Ycpr;
            validInd(1:Rx.CPR.Ntrain) = []; % discard symbols used in training carrier phase recovery algorithm
            
            % Phase tracking
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

       % Constellation plots
       if sim.Plots.isKey('Constellations') && sim.Plots('Constellations')
           figure(203), clf
           subplot(121)
           plot_constellation(Yeq(1, sim.Ndiscard+1:end-sim.Ndiscard), dataTX(1, Rx.AdEq.Ntrain+sim.Ndiscard+1:end-sim.Ndiscard), sim.M);
           axis square
           title('Pol X: After equalization')
           subplot(122)
           plot_constellation(Ycpr(1, sim.Ndiscard+1:end-sim.Ndiscard), dataTX(1, Rx.AdEq.Ntrain+Rx.CPR.Ntrain+sim.Ndiscard+1:end-sim.Ndiscard), sim.M);
           axis square
           title('Pol X: After CPR')   
           drawnow
       end
        
    elseif strcmpi(sim.ModFormat, 'DPSK')
%         [Y, Foff_est] = qpsk_freq_rec(Y, Rx.FreqRec, sim);
%         Y(:, 1:Rx.FreqRec.Ntrain) = []; % discard symbols used in training phase tracking algorithm     
%         validInd(1:Rx.FreqRec.Ntrain) = []; 
          Y = Yeq;
    end
    
    if strcmpi(sim.ModFormat, 'QAM')
        % Determines initial phase rotation for differential encoding
        dataRX = [demodulate(Y(1, :)); demodulate(Y(2, :))];
    elseif strcmpi(sim.ModFormat, 'DPSK')      
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
    ber.theory(k) = berAWGN(SNRdB_theory(k))
    1;
end

plots

if isfield(sim, 'Plots') && sim.Plots('BER') && length(ber.count) > 1
    figure(100), box on, hold on
    [~, link_attdB] = Fiber.link_attenuation(Tx.Laser.lambda);
    Prx = Tx.PlaunchdBm - link_attdB;
    plot(Prx, log10(ber.theory), '-')
    plot(Prx, log10(ber.count), '-o')
    legend('Theory including noise enhancement', 'Counted')
    xlabel('Received Power (dBm)')
    ylabel('log_{10}(BER)')
    axis([Prx(1) Prx(end) -8 0])
end
    
