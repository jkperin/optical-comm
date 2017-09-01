function [Lopt, Signal] = optimize_edf_length(E, Pump, Signal, Pon, SpanAttdB, verbose)
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
[Lopt, ~, exitflag] = fminbnd(@(L) -sum(objective(L, E, Pump, Signal, Pon, SpanAttdB)), 0, E.maxL);

if exitflag ~= 1
    warning('optimize_edf_length: optmization exited with exitflag %d\n', exitflag)
end

% Evaluate objective at optimal EDF length
onChs = objective(Lopt, E, Pump, Signal, Pon, SpanAttdB);
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
    

function onChs = objective(L, E, Pump, Signal, Pon, SpanAttdB)
    %% For a given EDF length L, determine maximum number of channels that can be on while satisfying GaindB >= SpanAttdB
    % Definitions
    maxIteration = 100;
    Poff = eps;
    
    % EDF fiber length
    E.L = L;
    
    % Start by setting all channels to a small power level
    Signal.P = Poff;
        
    GaindB = E.semi_analytical_gain(Pump, Signal);
    if all(GaindB < SpanAttdB) % no channel had gain higher than span attenuation
        onChs = 0;
        return;
    else % at least one channel has gain higher than span attenuation
        % Turn on channels that had gain above span attenuation
        onChs = (GaindB >= SpanAttdB);
        Signal.P(onChs) = Pon;
        n = 1;
        while n < maxIteration
            GaindB = E.semi_analytical_gain(Pump, Signal);
            onChsNew = (GaindB >= SpanAttdB);
            if not(isempty(onChsNew)) && all(onChs == onChsNew) % nothing changed
                return
            elseif not(isempty(onChsNew)) && all((onChs - onChsNew) < 0) 
                % 
                return 
            else
                % Turn off channel with smallest gain
                minGaindB = min(GaindB(onChs));
                Signal.P(GaindB == minGaindB) = Poff; % not a problem if more than one is turned off
                onChs(GaindB == minGaindB) = 0;
            end
        end
        
        if n == maxIteration
            warning('optimize_edf_length: optimization exceed maximum number of iterations')
        end       
    end
end
end
