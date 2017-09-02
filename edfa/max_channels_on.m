function onChs = max_channels_on(L, E, Pump, Signal, Pon, SpanAttdB)
    %% For a given EDF length L, determine maximum number of channels that can be on while satisfying GaindB >= SpanAttdB
    % Definitions
    maxIteration = 100;
    Poff = eps;
    
    % EDF fiber length
    E.L = L;
    
    % Start by setting all channels to a small power level
    Signal.P(1:end) = Poff;
        
    GaindB = E.semi_analytical_gain(Pump, Signal);
    if all(GaindB < SpanAttdB) % no channel had gain higher than span attenuation
        onChs = logical(zeros(size(Signal.wavelength)));
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
            elseif not(isempty(onChsNew)) && all((onChs - onChsNew) <= 0) 
                % 
                return 
            else
                % Turn off channel with smallest gain
                minGaindB = min(GaindB(onChs));
                Signal.P(GaindB == minGaindB) = Poff; % not a problem if more than one is turned off
                onChs(GaindB == minGaindB) = 0;
            end
            n = n + 1;
        end
        
        if n == maxIteration
            warning('optimize_edf_length: optimization exceed maximum number of iterations')
        end       
    end
end
