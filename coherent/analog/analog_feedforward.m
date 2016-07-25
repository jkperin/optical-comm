function [Xs, Analog] = analog_feedforward(Ys, Analog, sim, verbose)
%% Analog feedforward
% Inputs:
% - Ys: received signal after filtering to remove noise
% - totalLineWidth: sum of linewidths of transmitter and LO lasers
% - Analog: structure containing parameters of analog components, which
% include Mixer, Adder, Comparator
% - sim: simulation parameters {t: time, f: frequency, fs: sampling
% frequency}
% Outputs:
% - Xs: signal after four quadrant multiplier. Group delay is removed so
% sampling can be done without any offset
% - S: loop filter input
% - Sf: loop filter output
% - Analog: loop filter parameters

% Create components
% Mixers and adders for downconversion of pols X and Y
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

% Mixer and cubing for phase estimation
MixerIdQx = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
MixerQdIx = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
MixerIdQy = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
MixerQdIy = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
AdderX = AnalogAdder(Analog.Adder.filt, Analog.Adder.N0, sim.fs);
AdderY = AnalogAdder(Analog.Adder.filt, Analog.Adder.N0, sim.fs);

Cube1 = AnalogCubing(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
Cube2 = AnalogCubing(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
Cube3 = AnalogCubing(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
Cube4 = AnalogCubing(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);

% Low-pass filter
Fx = ClassFilter(Analog.Feedforward.LPF.filt, sim.fs);
Fy = ClassFilter(Analog.Feedforward.LPF.filt, sim.fs);
Fx.removeGroupDelay = false;
Fy.removeGroupDelay = false;

% Mixers and filters for regenerative frequency dividers
FreqDivXM1 = AnalogMixer(Analog.Feedforward.FreqDiv1.Mixer.filt, Analog.Feedforward.FreqDiv1.Mixer.N0, sim.fs);
FreqDivYM1 = AnalogMixer(Analog.Feedforward.FreqDiv1.Mixer.filt, Analog.Feedforward.FreqDiv2.Mixer.N0, sim.fs);
FreqDivXM2 = AnalogMixer(Analog.Feedforward.FreqDiv2.Mixer.filt, Analog.Feedforward.FreqDiv1.Mixer.N0, sim.fs);
FreqDivYM2 = AnalogMixer(Analog.Feedforward.FreqDiv2.Mixer.filt, Analog.Feedforward.FreqDiv2.Mixer.N0, sim.fs);

FreqDivXF1 = ClassFilter(Analog.Feedforward.FreqDiv1.filt, sim.fs);
FreqDivYF1 = ClassFilter(Analog.Feedforward.FreqDiv1.filt, sim.fs);
FreqDivXF2 = ClassFilter(Analog.Feedforward.FreqDiv2.filt, sim.fs);
FreqDivYF2 = ClassFilter(Analog.Feedforward.FreqDiv2.filt, sim.fs);

% Set filter memory to one so frequency divider can start
FreqDivXF1.memForward = ones(size(FreqDivXF1.memForward));
FreqDivYF1.memForward = ones(size(FreqDivYF1.memForward));
FreqDivXF2.memForward = ones(size(FreqDivXF2.memForward));
FreqDivYF2.memForward = ones(size(FreqDivYF2.memForward));

% Initializations
Y = [real(Ys(1, :));
     imag(Ys(1, :));
     real(Ys(2, :));
     imag(Ys(2, :))];

% Delay
freqDivDelay = max(Analog.Feedforward.freqDivDelay, 1);
Delay = round((Cube1.groupDelay + MixerIdQx.groupDelay + AdderX.groupDelay... % phase estimation
    + FreqDivXM1.groupDelay + FreqDivXF1.groupDelay... % regenerative frequency divider 1
    + FreqDivXM2.groupDelay + FreqDivXF2.groupDelay)*sim.fs); % regenerative frequency divider 2

% 4-th power phase estimation
% = Im{(xi + 1jxq)^4} = 4*(xi^3xq - xq^3xi)
Sx = 1/4*AdderX.add(MixerIdQx.mix(Cube1.cube(Y(1, :)), Y(2, :)), -MixerQdIx.mix(Cube2.cube(Y(2, :)), Y(1, :)));
Sy = 1/4*AdderY.add(MixerIdQy.mix(Cube3.cube(Y(3, :)), Y(4, :)), -MixerQdIy.mix(Cube4.cube(Y(4, :)), Y(3, :)));

% Low pass filtering and -1 gain
Sxf = -Fx.filter(Sx);
Syf = -Fy.filter(Sy);

% Initialization
Sxf = cos(2*pi*6e9*sim.t);
V2sinx = ones(size(Sxf));
V2siny = Syf;
Vsinx = V2sinx;
Vsiny = V2siny;
for t = freqDivDelay+1:length(sim.t)    
    % 1st regenerative frequency divider
%     V2sinx(t) = FreqDivXF1.filter(FreqDivXM1.mix(Sxf(t), V2sinx(t-freqDivDelay)));
    V2sinx(t) = 2*FreqDivXF1.filter(Sxf(t)*V2sinx(t-freqDivDelay));
    V2siny(t) = FreqDivYF1.filter(FreqDivYM1.mix(Syf(t), V2siny(t-freqDivDelay)));
    
    % 2nd regenerative frequency divider
    Vsinx(t) = FreqDivXF2.filter(FreqDivXM2.mix(V2sinx(t), Vsinx(t-freqDivDelay)));
    Vsiny(t) = FreqDivYF2.filter(FreqDivYM2.mix(V2siny(t), Vsiny(t-freqDivDelay)));    
end

% 90deg phase shift
Hshift = ifftshift(exp(-1j*pi/2*sign(sim.f)));
Vcosx = real(ifft(fft(Vsinx).*Hshift));
Vcosy = real(ifft(fft(Vsinx).*Hshift));

% Delay inputs
Y = circshift(Y, [0, -Delay]);

% Downconversion                 
X(1, :) = Sx1.add(Mx1.mix(Y(1, :), Vcosx), Mx2.mix(Y(2, :), Vsinx)); % pol X, I
X(2, :) = Sx2.add(-Mx3.mix(Y(1, :), Vsinx), Mx4.mix(Y(2, :), Vcosx)); % pol X, Q
X(3, :) = Sy1.add(My1.mix(Y(3, :), Vcosy), My2.mix(Y(4, :), Vsiny)); % pol Y, I
X(4, :) = Sy2.add(-My3.mix(Y(3, :), Vsiny), My4.mix(Y(4, :), Vcosy)); % pol Y, Q

% Remove group delay from signal path
delayDownconversion = sim.fs*(Mx1.groupDelay + Sx1.groupDelay); % Group delay in signal path
Hdelay = ifftshift(exp(1j*2*pi*sim.f/sim.fs*(Delay + delayDownconversion)));
X(1, :) = real(ifft(fft(X(1, :)).*Hdelay));
X(2, :) = real(ifft(fft(X(2, :)).*Hdelay));
X(3, :) = real(ifft(fft(X(3, :)).*Hdelay));
X(4, :) = real(ifft(fft(X(4, :)).*Hdelay));

% Build output
Xs = [X(1, :) + 1j*X(2, :); X(3, :) + 1j*X(4, :)];

if exist('verbose', 'var') && verbose
    figure(404), clf
%     subplot(212), hold on
    plot(Sx, 'DisplayName', 'pol X: phase estimation (PE)')
    plot(Sy, 'DisplayName', 'pol Y: PE')
    plot(Sx, 'DisplayName', 'pol X: PE fitlered')
    plot(Sy, 'DisplayName', 'pol Y: PE filtered')
    legend('-DynamicLegend')
end


