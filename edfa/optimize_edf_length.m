function [Lopt, Signal] = optimize_edf_length(E, Pump, Signal, Pon, SpanAttdB, method, verbose)
%% Optimize EDF length to maximize number of channels for which the EDF gain is larger than the span attenuation
% Inputs
% - E: isntance of class EDF
% - Pump: instance of class Channels corresponding to pump
% - Signal: instance of class Channels corresponding to signals
% - Pon: power of on channels
% - SpanAttdB: span attenuation at all signal wavelengths in dB
% Outputs:
% - Lopt: optimal EDF length
% - Signal: instance of class Channels corresponding to signals. The 

% Optimization
if strcmpi(method, 'fminbnd') % find optimal using Matlab's fminbnd function (Doesn't work consistently)
    [Lopt, ~, exitflag] = fminbnd(@(L) -sum(max_channels_on(L, E, Pump, Signal, Pon, SpanAttdB)), 0, E.maxL);

    if exitflag ~= 1
        warning('optimize_edf_length: optmization exited with exitflag %d\n', exitflag)
    end
elseif strcmpi(method, 'interp') % find optimal by interpolation
    Nsteps = 40;
    Ledf = linspace(1, E.maxL, Nsteps); % EDF is at least a meter long
    BW = zeros(size(Ledf));
    parfor k = 1:Nsteps
        BW(k) = sum(max_channels_on(Ledf(k), E, Pump, Signal, Pon, SpanAttdB));
    end
    
    Lfit = linspace(1, E.maxL, 1e3);
    [~, idx] = max(spline(Ledf, BW, Lfit));
    Lopt = Lfit(idx);
    
    if exist('verbose', 'var') && verbose
        figure(202), clf, hold on, box on
        stem(Ledf, BW)
        plot(Lopt, max(BW), 'or')
        xlabel('EDF length (m)')
        ylabel('Numer of channels On')
    end
end

% Evaluate objective at optimal EDF length
onChs = max_channels_on(Lopt, E, Pump, Signal, Pon, SpanAttdB);
Signal.P = 0;
Signal.P(onChs) = Pon;

% Plot
if exist('verbose', 'var') && verbose
    % Calculate optimal gain
    E.L = Lopt;
    Signal.P(not(onChs)) = eps; % small power just to measure gain
    GdB = E.semi_analytical_gain(Pump, Signal);    
    Signal.P(not(onChs)) = 0; % turn off again

    if length(SpanAttdB) == 1
        SpanAttdB = SpanAttdB*ones(size(Signal.wavelength));
    end
    
    figure(203), clf, hold on, box on
    plot(Signal.wavelength*1e9, GdB, 'DisplayName', 'EDF gain')
    plot(Signal.wavelength(onChs)*1e9, GdB(onChs), 'ob', 'DisplayName', 'Channels on')
    plot(Signal.wavelength(not(onChs))*1e9, GdB(not(onChs)), 'xr', 'DisplayName', 'Channels off')
    plot(Signal.wavelength*1e9, SpanAttdB, 'k', 'DisplayName', 'Span attenuation')
    xlabel('Wavelength (nm)')
    ylabel('Gain (dB)')
    legend('-DynamicLegend', 'Location', 'SouthEast')
    title(sprintf('EDF %s, L = %.2f m', E.type, E.L), 'Interpreter', 'none')
end
    

end
