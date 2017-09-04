function [E, Signal] = optimize_power_load_and_edf_length(method, E, Pump, Signal, problem, verbose)
%% Optimize power loading and EDF length
% Inputs
% - method: optimization method. 
% Assuming that channels can be either on or off: method = 'interp' or 'fminbnd'
% Assuming that channels power is upper bounded by Pon: method = 'particle swarm', or 'genetic'
% - E: isntance of class EDF
% - Pump: instance of class Channels corresponding to pump
% - Signal: instance of class Channels corresponding to signals
% - problem: struct with problem parameters: Pon, Namp, spanAttdB, and df
% Outputs:
% - Lopt: optimal EDF length
% - Signal: instance of class Channels corresponding to signals. The 

% Unpack problem parameters
Pon = problem.Pon;
Namp = problem.Namp;
spanAttdB = problem.spanAttdB;
df = problem.df;

%% Optimization
ploading = false;
switch lower(method)
    case 'fminbnd'
        %% Find EDF length that maximizes the number of ON channels
        % This uses Matlab's fminbnd function to search for the optimal fiber length
        % (Doesn't produce very consisten results)
        
        [Lopt, ~, exitflag] = fminbnd(@(L) -sum(max_channels_on(L, E, Pump, Signal, Pon, spanAttdB)), 0, E.maxL);

        if exitflag ~= 1
            warning('optimize_edf_length: optmization exited with exitflag %d\n', exitflag)
        end

        E.L = Lopt;

        % Evaluate objective at optimal EDF length
        onChs = max_channels_on(Lopt, E, Pump, Signal, Pon, spanAttdB);
        Signal.P = 0;
        Signal.P(onChs) = Pon;

    case 'interp'
        %% Find EDF length that maximizes the number of ON channels
        % Optimal EDF length is found by interpolation. Number of ON
        % channels is calculated for "Nsteps" points from EDF.L = 1 to
        % EDF.Lmax. Optimal result is obtained by interpolatin
        Nsteps = 40;
        Ledf = linspace(1, E.maxL, Nsteps); % EDF is at least a meter long
        BW = zeros(size(Ledf));
        parfor k = 1:Nsteps
            BW(k) = sum(max_channels_on(Ledf(k), E, Pump, Signal, Pon, spanAttdB));
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
        
        E.L = Lopt;
        
        % Evaluate objective at optimal EDF length
        onChs = max_channels_on(Lopt, E, Pump, Signal, Pon, spanAttdB);
        Signal.P = 0;
        Signal.P(onChs) = Pon;
        
    case 'particle swarm'
%         [E, Signal] = optimize_power_load_and_edf_length('interp', E, Pump, Signal, problem, false);
%         M = [E.L Signal.P;... % optimized by interp
%         	 E.L Pon*ones(1, Signal.N);... % Pon
%              E.L 1e-6*Pon*ones(1, Signal.N)]; %Poff
        
        options = optimoptions('particleswarm', 'Display', 'iter', 'UseParallel', true,...
            'MaxStallTime', 60, 'MaxStallIterations', 100, 'SwarmSize', min(200, 10*(Signal.N+1)));
        la = zeros(1, Signal.N+1); % lower bound
        lb = [E.maxL Pon*ones(1, Signal.N)]; % upper bound
        X = particleswarm(@(X) -capacity_linear_regime(X, E, Pump, Signal, spanAttdB, Namp, df), Signal.N+1, la, lb, options);
        
        Lopt = X(1);
        Signal.P = X(2:end);
        
        E.L = Lopt;
        GaindB = E.semi_analytical_gain(Pump, Signal);
        Signal.P(GaindB <= spanAttdB) = 0; % turn off channels that don't meet the gain requirement
        
        ploading = true;
    case 'genetic'
        options = optimoptions('ga', 'Display', 'iter', 'UseParallel', true, 'MaxStallTime', 60, 'MaxStallIterations', 100);
        la = zeros(1, Signal.N+1); % lower bound
        lb = [E.maxL Pon*ones(1, Signal.N)]; % upper bound
        X = ga(@(X) -capacity_linear_regime(X, E, Pump, Signal, spanAttdB, Namp, df), Signal.N+1, [], [], [], [], la, lb, [], options);
        
        Lopt = X(1);
        Signal.P = X(2:end);
        
        E.L = Lopt;
        GaindB = E.semi_analytical_gain(Pump, Signal);
        Signal.P(GaindB <= spanAttdB) = 0; % turn off channels that don't meet the gain requirement
        
        ploading = true;
    otherwise
        error('optimize_power_load_and_edf_length: invalid method')
end

% Plot
if exist('verbose', 'var') && verbose    
    fprintf('Optimization method = %s\n', method)
    fprintf('- Optimal EDF length = %.2f\n', E.L)
    fprintf('- Number of channels ON = %d\n', sum(Signal.P ~= 0))
    
    % Calculate optimal gain
    offChs = (Signal.P == 0);
    Signal.P(offChs) = eps; % small power just to measure gain
    GdB = E.semi_analytical_gain(Pump, Signal);    
    Signal.P(offChs) = 0; % turn off again

    if length(spanAttdB) == 1
        spanAttdB = spanAttdB*ones(size(Signal.wavelength));
    end
    
    figure(203), hold on, box on
    plot(Signal.wavelength*1e9, GdB, 'DisplayName', 'EDF gain')
    plot(Signal.wavelength(not(offChs))*1e9, GdB(not(offChs)), 'ob', 'DisplayName', 'Channels on')
    plot(Signal.wavelength(offChs)*1e9, GdB(offChs), 'xr', 'DisplayName', 'Channels off')
    plot(Signal.wavelength*1e9, spanAttdB, 'k', 'DisplayName', 'Span attenuation')
    xlabel('Wavelength (nm)')
    ylabel('Gain (dB)')
    legend('-DynamicLegend', 'Location', 'SouthEast')
    title(sprintf('EDF %s, L = %.2f m', E.type, E.L), 'Interpreter', 'none')
    
    if ploading
        figure(204), 
        subplot(211), box on
        plot(Signal.wavelength*1e9, Signal.PdBm)
        xlabel('Wavelength (nm)')
        ylabel('Signal power (dBm)')
        axis([Signal.wavelength([1 end])*1e9 10*log10(0.1*Pon/1e-3) 0])
        
        subplot(212), box on
        plot(Signal.wavelength*1e9, Signal.PdBm + GdB)
        xlabel('Wavelength (nm)')
        ylabel('Signal power after amplifier (dBm)')
    end
    drawnow
end