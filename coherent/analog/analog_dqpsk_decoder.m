function [dataRX, Yr] = analog_dqpsk_decoder(X, dataTX, Analog, sim)

Mxi = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
Mxq = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
Myi = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
Myq = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);

Analog.Filter.filt = []; %design_filter('bessel', 5, sim.Rs/(sim.fs/2));

Fxi = AnalogFilter(Analog.Filter.filt, sim.fs);
Fxq = AnalogFilter(Analog.Filter.filt, sim.fs);
Fyi = AnalogFilter(Analog.Filter.filt, sim.fs);
Fyq = AnalogFilter(Analog.Filter.filt, sim.fs);

% Full wave rectifiers (absolute operation)
abs1 = AnalogABS(Analog.ABS.filt, Analog.ABS.N0, sim.fs);
abs2 = AnalogABS(Analog.ABS.filt, Analog.ABS.N0, sim.fs);
abs3 = AnalogABS(Analog.ABS.filt, Analog.ABS.N0, sim.fs);
abs4 = AnalogABS(Analog.ABS.filt, Analog.ABS.N0, sim.fs);

% Comparators
CompX = AnalogComparator(Analog.Comparator.filt, Analog.Comparator.N0, sim.fs, Analog.Comparator.Vcc);
CompY = AnalogComparator(Analog.Comparator.filt, Analog.Comparator.N0, sim.fs, Analog.Comparator.Vcc);


% -90deg phase shift
phaseShift = zeros(size(sim.f));
phaseShift(sim.f > 0) = -pi/2;
phaseShift(sim.f < 0) = +pi/2;
phaseShift = fftshift(phaseShift);

Xxi = real(X(1, :));
Xxq = imag(X(1, :));
Xxq90 = real(ifft(fft(Xxq).*exp(1j*phaseShift))); % -90deg phase shift
Xyi = real(X(2, :));
Xyq = imag(X(2, :));
Xyq90 = real(ifft(fft(Xyq).*exp(1j*phaseShift))); % -90deg phase shift

Yr = zeros(4, length(X));
ctrl = zeros(2, length(X));
for t = sim.Mct+1:length(X)
    Yr(1, t) = Fxi.filter(Mxi.mix(Xxi(t), Xxi(t-sim.Mct)));
    Yr(2, t) = Fxq.filter(Mxq.mix(Xxq90(t), Xxq(t-sim.Mct)));
    ctrl(1, t) = CompX.compare(abs1.abs(Yr(1, t)), abs2.abs(Yr(2, t)));
    
    Yr(3, t) = Fyi.filter(Myi.mix(Xyi(t), Xyi(t-sim.Mct)));
    Yr(4, t) = Fyq.filter(Myq.mix(Xyq90(t), Xyq(t-sim.Mct)));
    ctrl(2, t) = CompY.compare(abs3.abs(Yr(3, t)), abs4.abs(Yr(4, t)));
end

% Remove group delay
delay = round((Mxi.groupDelay + Fxi.groupDelay)*sim.fs); 
ctrl = [ctrl(:, delay+1:end) ctrl(:, 1:delay)]; % remove group delay
Yr = [Yr(:, delay+1:end) Yr(:, 1:delay)]; % remove group delay

% Sampling
ctrl = ctrl(:, sim.Mct/2:sim.Mct:end);
Yr = Yr(:, sim.Mct/2:sim.Mct:end);

% Detection:
% 0 -> 0, 1 -> pi/2, 2 -> -pi/2, 3 -> pi
% a positive angle means counter-clockwise rotation
dataRX = zeros(2, length(Yr));
dataRX(1, ctrl(1, :) > 0 & Yr(1, :) > 0) = 0;
dataRX(1, ctrl(1, :) > 0 & Yr(1, :) < 0) = 3;
dataRX(1, ctrl(1, :) < 0 & Yr(2, :) > 0) = 1;
dataRX(1, ctrl(1, :) < 0 & Yr(2, :) < 0) = 2;


dataRX(2, ctrl(2, :) > 0 & Yr(3, :) > 0) = 0;
dataRX(2, ctrl(2, :) > 0 & Yr(3, :) < 0) = 3;
dataRX(2, ctrl(2, :) < 0 & Yr(4, :) > 0) = 1;
dataRX(2, ctrl(2, :) < 0 & Yr(4, :) < 0) = 2;
1;


