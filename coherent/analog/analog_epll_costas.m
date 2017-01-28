function [Xs, Analog, S, Sf] = analog_epll_costas(Ys, totalLineWidth, Analog, sim, verbose)
%% Analog electric phase locked loop via Costas loop
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
% Mixers and adders for single-sindeband mixer
[Mx1, Mx2, Mx3, Mx4, My1, My2, My3, My4] = Analog.Mixer.copy();
[Sx1, Sx2, Sy1, Sy2] = Analog.Adder.copy();

% Comparators
[Comp1, Comp2, Comp3, Comp4] = Analog.Comparator.copy();

% Mixer and Adders for Costas phase estimator stage
[MIdQx, MQdIx, MIdQy, MQdIy] = Analog.Mixer.copy();
[AdderX, AdderY, AdderXY]  = Analog.Adder.copy();

% Converts delay to number of samples in order to avoid interpolation
additionalDelay = max(round(Analog.Delay*sim.fs), 1); % delay is at least one sample

% Calculate group delay
totalGroupDelay = Analog.Mixer.groupDelay + Analog.Adder.groupDelay... % Four quadrant multiplier
    + Analog.Comparator.groupDelay + Analog.Mixer.groupDelay + Analog.CPRNpol*Analog.Adder.groupDelay... % phase estimation    
    + additionalDelay/sim.fs; % Additional loop delay e.g., propagation delay (minimum is 1/sim.fs since simulation is done in discrete time)
fprintf('Total loop delay: %.3f ps (%.2f bits, %d samples)\n', totalGroupDelay*1e12, totalGroupDelay*sim.Rb, ceil(totalGroupDelay*sim.fs));
Analog.totalGroupDelay = totalGroupDelay;

% Optimize EPLL parameters
if not(isfield(Analog, 'wn')) % wn was not yet defined; calculate optimal wn
    Analog.wn = optimizePLL(Analog.csi, totalGroupDelay, totalLineWidth, Analog.CPRNpol, sim, sim.shouldPlot('Phase error variance'));
end
Analog.EPLL.nums = [2*Analog.csi*Analog.wn Analog.wn^2];
Analog.EPLL.dens = [1 0 0]; % descending powers of s
[Analog.EPLL.numz, Analog.EPLL.denz] = impinvar(Analog.EPLL.nums, Analog.EPLL.dens, sim.fs);
fprintf('Loop filter fn: %.3f GHz\n', Analog.wn/(2*pi*1e9));

% Loop filter
LoopFilter = ClassFilter(Analog.EPLL.numz, Analog.EPLL.denz, sim.fs);

%% Loop
X = zeros(4, length(sim.t));
Xd = zeros(4, length(sim.t));
S = zeros(size(sim.t));
Sf = zeros(size(sim.t));

Yxi = real(Ys(1, :));
Yxq = imag(Ys(1, :));
Yyi = real(Ys(2, :));
Yyq = imag(Ys(2, :));

for t = additionalDelay+1:length(sim.t)
    % VCO: generates VCO output
    Vcos = cos(Sf(t-additionalDelay));
    Vsin = sin(Sf(t-additionalDelay));

    % Downconversion             
    X(1, t) = Sx1.add(Mx1.mix(Yxi(t), Vcos), -Mx2.mix(Yxq(t), Vsin)); % pol X, I
    X(2, t) = Sx2.add(Mx3.mix(Yxi(t), Vsin), Mx4.mix(Yxq(t), Vcos)); % pol X, Q
    X(3, t) = Sy1.add(My1.mix(Yyi(t), Vcos), -My2.mix(Yyq(t), Vsin)); % pol Y, I
    X(4, t) = Sy2.add(My3.mix(Yyi(t), Vsin), My4.mix(Yyq(t), Vcos)); % pol Y, Q

    % Phase estimation
    Xd(1, t) = Comp1.compare(X(1, t), 0);
    Xd(2, t) = Comp2.compare(X(2, t), 0);
    Xd(3, t) = Comp3.compare(X(3, t), 0);
    Xd(4, t) = Comp4.compare(X(4, t), 0);

    Sx = AdderX.add(MIdQx.mix(Xd(1, t), X(2, t)), -MQdIx.mix(Xd(2, t), X(1, t)));
    Sy = AdderY.add(MIdQy.mix(Xd(3, t), X(4, t)), -MQdIy.mix(Xd(4, t), X(3, t)));

    if Analog.CPRNpol == 2
        S(t) = AdderXY.add(Sx, Sy)/2; % loop filter input
    else
        S(t) = Sx; % loop filter input
    end

    % Loop filter
    Sf(t) = LoopFilter.filter(S(t));  
end

% Remove group delay from signal path
delay = (Analog.Mixer.groupDelay + Analog.Adder.groupDelay); % Group delay in signal path
Hdelay = ifftshift(exp(1j*2*pi*sim.f*delay));
X(1, :) = real(ifft(fft(X(1, :)).*Hdelay ));
X(2, :) = real(ifft(fft(X(2, :)).*Hdelay ));
X(3, :) = real(ifft(fft(X(3, :)).*Hdelay ));
X(4, :) = real(ifft(fft(X(4, :)).*Hdelay ));

% Build output
Xs = [X(1, :) + 1j*X(2, :); X(3, :) + 1j*X(4, :)];

 % Phase error plot -------------------------------------------------------
if exist('verbose', 'var') && verbose
    figure(404), clf
    subplot(211), hold on, box on
    plot(sim.t, S)
    xlabel('time (s)')
    ylabel('Loop filter input')       
    subplot(212), hold on, box on
    plot(sim.t, Sf)
    p = polyfit(sim.t, Sf, 1);
    plot(sim.t, polyval(p, sim.t));
    legend('VCO phase', sprintf('Linear fit (for freq offset) = %.2f GHz ramp', p(1)/(2*pi*1e9)))
    xlabel('time (s)')
    ylabel('Phase (rad)')
end

