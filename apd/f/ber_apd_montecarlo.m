function [ber, mpam] = ber_apd_montecarlo(mpam, Tx, Fiber, Apd, Rx, sim)
%% Calculate BER of unamplified IM-DD system with APD detector using montecarlo simulation

% System received pulse shape frequency response
HrxPshape = apd_system_received_pulse_shape(mpam, Tx, Fiber, Apd, Rx, sim); % this is only used if rx.filtering = matched or eq.type = fixed...

% Ajust levels to desired transmitted power and extinction ratio
mpam = mpam.adjust_levels(Tx.Ptx, Tx.Mod.rexdB);
AGC = 1/(mpam.a(end)*Apd.Gain*Apd.R*Fiber.link_attenuation(Tx.Laser.wavelength));

%% Modulated PAM signal
% dataTX = [0 0 0 0 0 1 0 0 0 0 0];
dataTX = randi([0 mpam.M-1], 1, sim.Nsymb); % Random sequence
Nzero = 10;
dataTX([1:Nzero end-Nzero+1:end]) = 0;
xk = mpam.signal(dataTX);

%% DAC
xt = dac(xk, Tx.DAC, sim, sim.shouldPlot('DAC output')); 

%% Generate optical signal
[Etx, Pt] = eam(Tx.Laser.cw(sim), xt, Tx.Mod, sim.f);

%% Ensures that transmitted power is at the right level
AGC = AGC/(Tx.Ptx/mean(Pt));
Etx = Etx*sqrt(Tx.Ptx/mean(Pt));

%% Fiber propagation
Erx = Fiber.linear_propagation(Etx, sim.f, Tx.Laser.wavelength);

%% Detect and add noises
% yt = Apd.detect(Erx, sim.fs, 'no noise');
yt = Apd.detect(Erx, sim.fs, 'gaussian', Rx.N0);

%% Whitening filter
if sim.WhiteningFilter
    [~, yt] = Apd.Hwhitening(sim.f, mean(abs(Erx).^2), Rx.N0, yt);
end

%% Automatic gain control
% Normalize signal so that highest level is equal to 1
mpam = mpam.norm_levels;
yt = yt*AGC;
yt = yt - mean(yt) + mean(mpam.a);

%% ADC
% ADC performs filtering, quantization, and downsampling
% For an ideal ADC, ADC.ENOB = Inf
% Align received and transmitted signals
% Rx.ADC.offset = 0;
Rx.ADC.timeRefSignal = xt; % align filtered signal ytf to this reference
switch lower(Rx.filtering)
    case 'antialiasing' % receiver filter is specified in ADC.filt
        [yk, ~, ytf] = adc(yt, Rx.ADC, sim);
    case 'matched' % receiver filter is matched filter
        Hrx = conj(HrxPshape);
        [yk, ~, ytf] = adc(yt, Rx.ADC, sim, Hrx); % replace ADC antialiasing filter by matched filter
    otherwise
        error('ber_preamp_sys_montecarlo: Rx.filtering must be either antialiasing or matched')
end     

%% Equalization
Rx.eq.trainSeq = dataTX;
[yd, Rx.eq] = equalize(Rx.eq, yk, HrxPshape, mpam, sim);

% Symbols to be discard in BER calculation
ndiscard = [1:Rx.eq.Ndiscard(1)+sim.Ndiscard (sim.Nsymb-Rx.eq.Ndiscard(2)-sim.Ndiscard):sim.Nsymb];
ydfull = yd;
yd(ndiscard) = []; 
dataTX(ndiscard) = [];

%% Demodulate
if mpam.optimize_level_spacing
    dataRX = mpam.demod(yd);
else
    [dataRX, mpam] = mpam.demod_sweeping_thresholds(yd, dataTX);
end

%% True BER
[~, ber] = biterr(dataRX, dataTX);

%% Plots
if sim.shouldPlot('Empirical noise pdf')
    % Empirical pdf for a level
    figure(100)
    [nn, xx] = hist(yd(dataTX == 2), 50);
    nn = nn/trapz(xx, nn);
    bar(xx, nn)
    title('Empirical pdf for PAM level 2')
end    

if sim.shouldPlot('Equalizer')
    figure(101), clf
    for k = 1:size(Rx.eq.h, 2)
        [h, w] = freqz(Rx.eq.h(:, k), 1);
        subplot(121), hold on, box on
        plot(w/(2*pi), abs(h).^2)
        
        subplot(122), hold on, box on
        plot(w/(2*pi), unwrap(angle(h)))
    end
    subplot(121), hold on, box on
    xlabel('Normalized frequency')
    ylabel('|H(f)|^2')
    title('Equalizer amplitude response')

    subplot(122), hold on, box on
    xlabel('Normalized frequency')
    ylabel('arg(H(f))')
    title('Equalizer phase response')
    drawnow
end

if sim.shouldPlot('Optical eye diagram')  
    Ntraces = 500;
    Nstart = sim.Ndiscard*sim.Mct + 1;
    Nend = min(Nstart + Ntraces*2*sim.Mct, length(Etx));
    figure(103), clf, box on
    eyediagram(abs(Etx(Nstart:Nend)).^2, 2*sim.Mct)
    title('Optical eye diagram')
    drawnow
end

if sim.shouldPlot('Received signal eye diagram')
    mpam = mpam.norm_levels();
    Ntraces = 500;
    Nstart = sim.Ndiscard*sim.Mct + 1;
    Nend = min(Nstart + Ntraces*2*sim.Mct, length(ytf));
    figure(104), clf, box on, hold on
    eyediagram(ytf(Nstart:Nend), 2*sim.Mct)
    title('Received signal eye diagram')
    a = axis;
    h1 = plot(a(1:2), mpam.a*[1 1], '-k');
    h2 = plot(a(1:2), mpam.b*[1 1], '--k');
    h3 = plot((sim.Mct+1)*[1 1], a(3:4), 'k');
    legend([h1(1) h2(1) h3], {'Levels', 'Decision thresholds', 'Sampling point'})
    drawnow
end

if sim.shouldPlot('Signal after equalization')
    mpam = mpam.norm_levels();
    figure(105), clf, box on, hold on
    h1 = plot(ydfull, 'o');
    a = axis;
    h2= plot(a(1:2), mpam.a*[1 1], '-k');
    h3 = plot(a(1:2), mpam.b*[1 1], '--k');
    h4 = plot(Rx.eq.Ndiscard(1)*[1 1], a(3:4), ':k');
    h5 = plot((sim.Nsymb-Rx.eq.Ndiscard(2))*[1 1], a(3:4), ':k');
    legend([h1 h2(1) h3(1) h4], {'Equalized samples', 'PAM levels',...
        'Decision thresholds', 'BER measurement window'})
    title('Signal after equalization')
    axis([1 sim.Nsymb -0.2 1.2])
    drawnow
end
