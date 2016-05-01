function [ber, SNRdBtheory] = ber_coherent_analog_epll(Tx, Fiber, Rx, sim)
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
berAWGN = @(SNRdB) berawgn(SNRdB - 3 - 10*log10(log2(M)), lower(sim.ModFormat), M);
% Note: -3 is due to the fact SNR = 2Es/N0

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
noiseBW = sim.Rs/2; %trapz(sim.f, abs(eq.Hrx.*eq.Hff(sim.f/sim.Rs)).^2)/2;
sim.f = sim.f.';

Analog = Rx.Analog;

% Create components
% Mixers and adders for four quadrant multiplier for X and Y
Mx1 = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
Mx2 = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
Mx3 = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
Mx4 = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
Sx1 = AnalogAdder(Analog.Adder.filt, Analog.Adder.N0, sim.fs);
Sx2 = AnalogAdder(Analog.Adder.filt, Analog.Adder.N0, sim.fs);

My1 = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
My2 = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
My3 = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
My4 = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
Sy1 = AnalogAdder(Analog.Adder.filt, Analog.Adder.N0, sim.fs);
Sy2 = AnalogAdder(Analog.Adder.filt, Analog.Adder.N0, sim.fs);

% Comparators
Comp1 = AnalogComparator(Analog.Comparator.filt, Analog.Comparator.N0, sim.fs, Analog.Comparator.Vcc);
Comp2 = AnalogComparator(Analog.Comparator.filt, Analog.Comparator.N0, sim.fs, Analog.Comparator.Vcc);
Comp3 = AnalogComparator(Analog.Comparator.filt, Analog.Comparator.N0, sim.fs, Analog.Comparator.Vcc);
Comp4 = AnalogComparator(Analog.Comparator.filt, Analog.Comparator.N0, sim.fs, Analog.Comparator.Vcc);

MixerIdQx = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
MixerQdIx = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
MixerIdQy = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
MixerQdIy = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);

AdderX = AnalogAdder(Analog.Adder.filt, Analog.Adder.N0, sim.fs);
AdderY = AnalogAdder(Analog.Adder.filt, Analog.Adder.N0, sim.fs);
AdderXY = AnalogAdder(Analog.Adder.filt, Analog.Adder.N0, sim.fs);

% Calculate group delay
totalGroupDelay = Mx1.groupDelay + Sx1.groupDelay... % Four quadrant multiplier
    + Comp1.groupDelay + MixerIdQx.groupDelay + AdderX.groupDelay + AdderXY.groupDelay... % phase estimation    
    + Analog.Delay/sim.fs; % Additional loop delay e.g., propagation delay (minimum is 1/sim.fs since simulation is done in discrete time)
fprintf('Total loop delay: %.3f ps (%.2f bits, %d samples)\n', totalGroupDelay*1e12, totalGroupDelay*sim.Rb, ceil(totalGroupDelay*sim.fs));

% Optimize EPLL parameters
Analog.wn = optimizePLL(Analog.csi, Analog.Kdc, totalGroupDelay, Tx.Laser.linewidth+Rx.LO.linewidth, sim);
Analog.EPLL.nums = Analog.Kdc*[2*Analog.csi*Analog.wn Analog.wn^2];
Analog.EPLL.dens = [1 0 0]; % descending powers of s
[Analog.EPLL.numz, Analog.EPLL.denz] = impinvar(Analog.EPLL.nums, Analog.EPLL.dens, sim.fs);
fprintf('Loop filter fn: %.3f GHz\n', Analog.wn/(2*pi*1e9));

%
numz = Analog.EPLL.numz;
denz = Analog.EPLL.denz;
numLen = length(numz);
denLen = length(denz);

% Transmitted power swipe
ber.count = zeros(size(Tx.PlaunchdBm));
ber.theory = zeros(size(Tx.PlaunchdBm));
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

    % Ys = [Yxi + 1j*Yxq; Yyi + 1j*Yyq];
          
    %% Pol demux
    %% !!! Must be implemented 
    
    % Estimate SNR including noise enhacement penalty of an ideal symbol-rate linear equalizer   
    Prx = dBm2Watt(Tx.Laser.PdBm)/Fiber.link_attenuation(Tx.Laser.lambda); % received power
    Plo = dBm2Watt(Rx.LO.PdBm);                                            % local oscillator power
    Ppd = abs(sqrt(Plo/(4*sim.Npol)) + sqrt(Prx/(4*sim.Npol))).^2;         % incident power in each photodiode
    Psig = Plo*Prx/(2*sim.Npol*sim.Npol);                          % Signal power per real dimension
    varShot = 2*Rx.PD.varShot(Ppd, noiseBW);                               % Shot noise variance per real dimension
    varThermal = Rx.N0*noiseBW;                                            % Thermal noise variance per real dimension
    SNRdBtheory(k) = 10*log10(2*Psig/(varShot + varThermal));

    %% Carrier phase recovery algorithm  
    X = zeros(4, length(sim.t));
    Xd = zeros(4, length(sim.t));
    S = zeros(size(sim.t));
    Sf = zeros(size(sim.t));
       
    if strcmpi(sim.ModFormat, 'QAM') % not done if DQPSK
        for t = max([Analog.Delay sim.Mct numLen denLen])+1:length(sim.t)
            % VCO: generates VCO output
            Vcos = cos(Sf(t-Analog.Delay));
            Vsin = sin(Sf(t-Analog.Delay));
            
            % Four quadrant multiplier                  
            X(1, t) = Sx1.add(Mx1.mix(Yxi(t), Vcos), -Mx2.mix(Yxq(t), Vsin)); % pol X, I
            X(2, t) = Sx2.add(Mx3.mix(Yxi(t), Vsin), Mx4.mix(Yxq(t), Vcos)); % pol X, Q
            X(3, t) = Sy1.add(My1.mix(Yyi(t), Vcos), -My2.mix(Yyq(t), Vsin)); % pol Y, I
            X(4, t) = Sy2.add(My3.mix(Yyi(t), Vsin), My4.mix(Yyq(t), Vcos)); % pol Y, Q

            % Phase estimation
            Xd(1, t) = Comp1.compare(X(1, t), Analog.Comparator.Vref);
            Xd(2, t) = Comp2.compare(X(2, t), Analog.Comparator.Vref);
            Xd(3, t) = Comp3.compare(X(3, t), Analog.Comparator.Vref);
            Xd(4, t) = Comp4.compare(X(4, t), Analog.Comparator.Vref);
                
            Sx = AdderX.add(MixerIdQx.mix(Xd(1, t), X(2, t)), -MixerQdIx.mix(Xd(2, t), X(1, t)));
            Sy = AdderY.add(MixerIdQy.mix(Xd(3, t), X(4, t)), -MixerQdIy.mix(Xd(4, t), X(3, t)));
            S(t) = AdderXY.add(Sx, Sy);

            % Loop filter
            Sf(t) = sum(numz.*S(t:-1:t-numLen+1))...
                - sum(Sf(t-1:-1:t-denLen+1).*denz(2:end));  
        end
        
        % Reset components states
        Mx1.reset(); Mx2.reset(); Mx3.reset(); Mx4.reset();
        My1.reset(); My2.reset(); My3.reset(); My4.reset();
        Sx1.reset(); Sx2.reset();
        Sy1.reset(); Sy2.reset();
        Comp1.reset(); Comp2.reset(); Comp3.reset(); Comp4.reset();
        MixerIdQx.reset(); MixerQdIx.reset();
        MixerIdQy.reset(); MixerQdIy.reset();
        AdderX.reset();
        AdderY.reset();
        AdderXY.reset();
        
    elseif strcmpi(sim.ModFormat, 'DPSK')
        1;
    end
    
    % Build output
    Xs = sqrt(sqrt(2))*[X(1, :) + 1j*X(2, :); X(3, :) + 1j*X(4, :)];
  
    % Sampling
    delay = round((Mx1.groupDelay + Sx1.groupDelay)*sim.fs); % Group delay of four quadrant multiplier
    Xs = circshift(Xs, -delay, 2); % remove group delay
    Y = Xs(:, sim.Mct/2:sim.Mct:end);  
    % Note: decisions could be obtained directly from Xd. Xs is used here 
    % instead so we can see the constellation 
    
    % Correct for constant phase rotation
    [Y, phiPT] = phase_estimation(Y, Rx.PT, symbolsTX(:, validInd), sim);
          
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
    ber.theory(k) = berAWGN(SNRdBtheory(k))

    % Phase error plot
    if sim.Plots.isKey('EPLL phase error') && sim.Plots('EPLL phase error')
        figure(204), clf
        subplot(211), hold on, box on
        plot(sim.t, S)
        xlabel('time (s)')
        ylabel('Loop filter input')       
        subplot(212), hold on, box on
        plot(sim.t, Sf)
        p = polyfit(sim.t, Sf, 1);
        plot(sim.t, polyval(p, sim.t));
        plot(sim.t, sign(p(1))*2*pi*sim.t*Rx.LO.freqOffset)
        legend('VCO phase','Linear fit (for freq offset)',  'Frequency offset ramp')
        xlabel('time (s)')
        ylabel('Phase (rad)')
    end
    
   % Constellation plots
   if strcmpi(sim.ModFormat, 'QAM') && sim.Plots.isKey('Constellations') && sim.Plots('Constellations')
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
    
