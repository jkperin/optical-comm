function ber = ber_coherent_dsp_matched_filtering(Tx, Fiber, Rx, sim)
%% Simulate transmission of coherent system. The receiver consists of 
%% matched filter matched to the channel impulse response, symbol-rate sampling
%% and linear equalization

ModFormat = sim.ModFormat;

% Ensures that modulation format is QPSK or DQPSK
assert(ModFormat.M == 4, 'ber_coherent_dsp: Modulation format must be quarternary')

%% Generate transmitted symbols
% Note 1: this implement assumes that no DAC is being used so ModFormat.signal
% generates rectangular pulses sampled at sim.fs
% Note 2: ModFormat.signal does not remove group delay due to pulse shaping
% so the center of the pulse is at (Mct+1)/2, if pulse shape is 'rect'
dataTX = randi([0 ModFormat.M-1], [2, sim.Nsymb]); % symbol stream for each polarization
Nzero = 10; % zero Nzero first and last symbols to make sequence periodic
dataTX(:, 1:Nzero) = 0;
dataTX(:, end-Nzero+1:end) = 0;

[Vin, symbolsTX] = ModFormat.signal(dataTX); % X & Y pol

% Filter drive waveforms for modulators txfilt.
% group delay of Tx.filt.H has already been removed
Htx = ifftshift(Tx.filt.H(sim.f/sim.fs).*exp(1j*2*pi*sim.f/sim.fs*ModFormat.pulse_shape_grpdelay)); % transmitter filter and remove group delay due to pulse shaping in ModFormat
Vout(1, :) = real(ifft(fft(real(Vin(1, :))).*Htx)) + 1j*real(ifft(fft(imag(Vin(1, :))).*Htx)); 
Vout(2, :)= real(ifft(fft(real(Vin(2, :))).*Htx)) + 1j*real(ifft(fft(imag(Vin(2, :))).*Htx));

% Define fixed time-domain symbol-rate linear equalizer 
% mpam is only used to get pulse shape
eq.type = 'fixed td-sr-le';
eq.Ntaps = 31;
HtxPshape = freqz(ones(1, sim.Mct)/sim.Mct, 1, sim.f, sim.fs).*exp(1j*2*pi*sim.f/sim.fs*(sim.Mct-1)/2);
HrxPshape = HtxPshape.*Tx.filt.H(sim.f/sim.fs).*Tx.Mod.H;

% BER in AWGN channel
ber.theory = ber_coherent_awgn(Tx, Fiber, Rx, sim);

%% Swipe launched power
% Transmitted power swipe
ber.count = zeros(size(Tx.PlaunchdBm));
validRange = sim.Ndiscard+1:sim.Nsymb-sim.Ndiscard;
M = ModFormat.M;
Enorm = mean(abs(pammod(0:log2(M)-1, log2(M))).^2);
for k = 1:length(Tx.PlaunchdBm)
    fprintf('-- Launch power: %.2f dBm\n', Tx.PlaunchdBm(k));   
    
    %% ========= Modulator ==========
    Tx.Laser.PdBm = Tx.PlaunchdBm(k);
    Ein = Tx.Laser.cw(sim); % Generates electric field with intensity and phase noise       
    Ein = mzm(Ein, Vout, Tx.Mod); % modulate optical signal using eletro-optical modulator (EOM)
    
    % Ensure that transmitted power is a desired level
    Ein = Ein*sqrt(dBm2Watt(Tx.PlaunchdBm(k))/sum(mean(abs(Ein).^2, 2)));

    %% ========= Propagation ========== 
    Erec = Fiber.linear_propagation(Ein, sim.f, Tx.Laser.lambda);
    
    %% ====== Optical Amplifier =======
    if isfield(sim, 'preAmp') && sim.preAmp % only included if sim.preAmp is true
        disp('- IMPORTANT: Simulation including optical amplifier!')
        [Erec, OSNRdBtheory] = Rx.OptAmp.amp(Erec, sim.fs);
        fprintf('OSNR = %.2f dB\n', OSNRdBtheory)
        
        % Adjust power to pre-defined value
        Att = dBm2Watt(Rx.OptAmpOutPowerdBm)/dBm2Watt(power_meter(Erec));
        Erec = Erec*sqrt(Att);  % keep constant received power
        
        % check
%         OSNRdBtheorycheck = 10*log10(dBm2Watt(Tx.PlaunchdBm(k))/(2*Rx.OptAmp.nsp*Rx.OptAmp.h*Rx.OptAmp.c/Tx.Laser.lambda*12.5e9))
%         SNRdBtheorycheck = 10*log10(dBm2Watt(Tx.PlaunchdBm(k))/(Rx.OptAmp.nsp*Rx.OptAmp.h*Rx.OptAmp.c/Tx.Laser.lambda*sim.Rs))
    end

    %% ========= Receiver =============
    ELO = Rx.LO.cw(sim); % generates continuous-wave electric field in 1 pol with intensity and phase noise
    ELO = [sqrt(1/2)*ELO;    % LO field in x polarization at PBS or BS input
           sqrt(1/2)*ELO];    % LO field in y polarization at PBS or BS input
    Y = dual_pol_coherent_receiver(Erec, ELO, Rx, sim);
        
    %% Matched filtering
    Hrx = ifftshift(conj(HrxPshape));
    Yxi = real(ifft(fft(real(Y(1, :))).*Hrx));
    Yxq = real(ifft(fft(imag(Y(1, :))).*Hrx));
    Yyi = real(ifft(fft(real(Y(2, :))).*Hrx));
    Yyq = real(ifft(fft(imag(Y(2, :))).*Hrx));
    
    %% Symbol-rate sampling
    Yxi = Yxi(1:sim.Mct:end);
    Yxq = Yxq(1:sim.Mct:end);
    Yyi = Yyi(1:sim.Mct:end);
    Yyq = Yyq(1:sim.Mct:end);
    
    %% Automatic gain control (AGC)
    Yxi = sqrt(Enorm./mean(abs(Yxi).^2))*Yxi;
    Yxq = sqrt(Enorm./mean(abs(Yxq).^2))*Yxq;
    Yyi = sqrt(Enorm./mean(abs(Yyi).^2))*Yyi;
    Yyq = sqrt(Enorm./mean(abs(Yyq).^2))*Yyq;
    
    %% Equalization
    [Ydxi, eq] = equalize(eq, Yxi, HrxPshape, [], sim);
    [Ydxq, eq] = equalize(eq, Yxq, HrxPshape, [], sim);
    [Ydyi, eq] = equalize(eq, Yyi, HrxPshape, [], sim);
    [Ydyq, eq] = equalize(eq, Yyq, HrxPshape, [], sim);    
        
    % Finer gain control: to make sure that QAM constellation is of the
    % right size
    Enorm = mean(abs(pammod(0:log2(ModFormat.M)-1, log2(ModFormat.M))).^2);
    Ydxi = sqrt(Enorm/mean(abs(Ydxi).^2))*Ydxi;
    Ydxq = sqrt(Enorm/mean(abs(Ydxq).^2))*Ydxq;
    Ydyi = sqrt(Enorm/mean(abs(Ydyi).^2))*Ydyi;
    Ydyq = sqrt(Enorm/mean(abs(Ydyq).^2))*Ydyq;
    
    Yd = [Ydxi + 1j*Ydxq; Ydyi + 1j*Ydyq];
        
    %% Demodulate
    Yd = Yd(:, validRange);
    if strcmpi(ModFormat.type, 'QAM')        
        dataRX = ModFormat.demod(Yd);
    elseif strcmpi(ModFormat.type, 'DPSK')
        dataRX = ModFormat.demod(Yd, symbolsTX(:, sim.Ndiscard));
    end
    
    %% BER calculation
    [~, berX(k)] = biterr(dataTX(1, validRange), dataRX(1, :))
    [~, berY(k)] = biterr(dataTX(2, validRange), dataRX(2, :))
    ber.count(k) = 0.5*(berX(k) + berY(k));

   %% Constellation plots
   if sim.shouldPlot('Constellations')
       figure(203), clf
       subplot(121)
       plot_constellation(Yd(1, :), dataTX(1, validRange), ModFormat.M);
       axis square
       subplot(122)
       plot_constellation(Yd(2, :), dataTX(2, validRange), ModFormat.M);
       axis square
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
    hline(1) = plot(Prx, log10(ber.theory), '-');
    hline(2) = plot(Prx, log10(ber.count), '-o');
    legend(hline, {'Theory', 'Counted'})
    xlabel('Received Power (dBm)')
    ylabel('log_{10}(BER)')
    axis([Prx(1) Prx(end) -8 0])
end
    
