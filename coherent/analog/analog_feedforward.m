function [Xs, Analog] = analog_feedforward(Ys, Analog, sim)
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

% Mixers and filters for regenerative frequency dividers
MdivX1 = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
MdivX2 = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
MdivY1 = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);
MdivY2 = AnalogMixer(Analog.Mixer.filt, Analog.Mixer.N0, sim.fs);

FdivX1 = AnalogFilter(Analog.FeedforwardLPF.filt, sim.fs);
FdivX2 = AnalogFilter(Analog.FeedforwardLPF.filt, sim.fs);
FdivY1 = AnalogFilter(Analog.FeedforwardLPF.filt, sim.fs);
FdivY2 = AnalogFilter(Analog.FeedforwardLPF.filt, sim.fs);

% -90deg phase shift filter
H = exp(-1j*pi/2*sign(sim.f));
h = real(fftshift(ifft(ifftshift(H))));
tol = 1e-3;
Eh = cumtrapz(sim.t, abs(h).^2)/trapz(sim.t, abs(h).^2);
n = -sim.N/2:sim.N/2-1;
n1 = find(Eh >= 1-tol, 1, 'first');
n2 = find(Eh <= tol, 1, 'last');
nmax = max(abs([n(n1) n(n2)]));
h = h(n > -nmax & n < nmax);
n = n(n > -nmax & n < nmax);
shifterLen = length(h);

filt.h = h;
filt.num = h;
filt.den = 1;
phaseShifterX =  AnalogFilter(filt, sim.fs);
phaseShifterY =  AnalogFilter(filt, sim.fs);

% Initializations
X = zeros(4, length(sim.t));

Sx = zeros(size(sim.t));
Sy = zeros(size(sim.t));

Yxi = real(Ys(1, :));
Yxq = imag(Ys(1, :));
Yyi = real(Ys(2, :));
Yyq = imag(Ys(2, :));

% Delay
freqDivDelay = round((MdivX1.groupDelay + FdivX1.groupDelay)*sim.fs);
shifterDelay = (shifterLen-1)/2; % delay of symmetric FIR filter
Delay = round((Cube1.groupDelay + MixerIdQx.groupDelay + AdderX.groupDelay... % phase estimation
    + MdivX1.groupDelay + FdivX1.groupDelay... % regenerative frequency divider 1
    + MdivX2.groupDelay + FdivX2.groupDelay)*sim.fs... % regenerative frequency divider 2
    + shifterDelay);

for t = 1:length(sim.t)
    % Phase estimation
    % = Im{(xi + 1jxq)^4} = 4*(xi^3xq - xq^3xi)
    Sx(t) = -1/4*AdderX.add(MixerIdQx.mix(Cube1.cube(Yxi(t)), Yxq(t)), -MixerQdIx.mix(Cube2.cube(Yxq(t)), Yxi(t)));
    Sy(t) = -1/4*AdderY.add(MixerIdQy.mix(Cube3.cube(Yyi(t)), Yyq(t)), -MixerQdIy.mix(Cube4.cube(Yyq(t)), Yyi(t)));
end

% Initialization
Vx = Sx;
Vsinx = Sx;
Vy = Sy;
Vsiny = Sy;
Vcosx = zeros(size(sim.t));
Vcosy = zeros(size(sim.t));
AA = cos(2*pi*1e9*sim.t);
Vx = ones(size(AA));
for t = max([Delay sim.Mct])+1:length(sim.t)
    Vx(t) = FdivX1.filter(2*MdivX1.mix(AA(t), Vx(t-freqDivDelay)));
    Vy(t) = FdivX2.filter(2*MdivY2.mix(Sy(t), Vy(t-freqDivDelay)));
    
    Vsinx(t) = regenerative_frequency_divider(Vx(t), Vsinx(t-freqDivDelay), MdivX2, FdivX2);
    Vsiny(t) = regenerative_frequency_divider(Vy(t), Vsiny(t-freqDivDelay), MdivY2, FdivY2);
    
    Vcosx(t) = -phaseShifterX.filter(Vsinx(t));
    Vcosy(t) = -phaseShifterY.filter(Vsiny(t));
    
    % Downconversion                 
    X(1, t) = Sx1.add(Mx1.mix(Yxi(t-Delay), Vcosx(t)), Mx2.mix(Yxq(t-Delay), Vsinx(t-shifterDelay))); % pol X, I
    X(2, t) = Sx2.add(-Mx3.mix(Yxi(t-Delay), Vsinx(t-shifterDelay)), Mx4.mix(Yxq(t-Delay), Vcosx(t))); % pol X, Q
    X(3, t) = Sy1.add(My1.mix(Yyi(t-Delay), Vcosy(t)), My2.mix(Yyq(t-Delay), Vsiny(t-shifterDelay))); % pol Y, I
    X(4, t) = Sy2.add(-My3.mix(Yyi(t-Delay), Vsiny(t-shifterDelay)), My4.mix(Yyq(t-Delay), Vcosy(t))); % pol Y, Q
end

% Build output
Xs = [X(1, :) + 1j*X(2, :); X(3, :) + 1j*X(4, :)];

% Sampling
delay = round((Cube1.grpdelay + MixerIdQx.grpdelay + AdderX.grpdelay... % phase estimation
    + MdivX1.groupDelay + FdivX1.grpdelay... % regenerative frequency divider 1
    + MdivX2.groupDelay + FdivX2.grpdelay)*sim.fs... % regenerative frequency divider 2
    + shifterDelay... % 90deg phase shift
    + (Mx1.groupDelay + Sx1.groupDelay)*sim.fs); % downconversion
Xs = [Xs(:, delay+1:end) Xs(:, 1:delay)]; % remove group delay

    function y = regenerative_frequency_divider(x, ypast, mix, filt)
        y = filt.filter(2*mix.mix(x, ypast));
    end
end