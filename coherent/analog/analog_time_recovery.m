function [ysamp, idx, Ndiscard] = analog_time_recovery(yct, TimeRec, sim, verbose)
%% Time recovery for analog-based receiver
% Inputs:
% - yc: input signal in "continuous-time". If yct is complex, ony the I 
% component will be used for estimating the clock.
% - TimeRec: struct containing time recovery method properties
%     - type: eithehr 'squarer' or 'Mueller-Muller'
%     - bpf (only if type = 'squarer'): band pass filter struct 
%     - csi (type = 'Mueller-Muller'): damping of loop filter
%     - wn (type = 'Mueller-Muller'): relaxation frequency of loop filter
%     - CT2DT (type = 'Mueller-Muller', optional, default = 'bilinear')
%     method of converting continuous-time loop filter into discrete time
% - sim: struct containing simulation parameters
% Outputs:
% - ysamp: sampled signal
% - idx: indices such that ysamp = yct(idx)

if not(isfield(TimeRec, 'Ndiscard'))
    Ndiscard = [0 0];
end

switch lower(TimeRec.type)
    case 'spectral-line-bpf'
        %% Spectral line syncronizer based on squaring and bandpass filtering
        % Differentiator -> nonlinearity (squaring) -> band pass filter
        % centered at the symbol rate
        % Note: BPF could be replaced by a PLL, but this is not considered
        % here
        
        yct = resample(yct, TimeRec.Mct, sim.Mct);
        fs = sim.fs*TimeRec.Mct/sim.Mct;
        [f, t] = freq_time(length(yct), fs);
        Squarer = AnalogSquaring(TimeRec.squarerFilt, TimeRec.N0, fs);
        
        % Only uses I component in estimating clock
        y = real(yct);
               
        % Nonlinearity
        ynl = Squarer.square(y);
        % Note: group delay due to filtering operations in squarer must be
        % removed completely to preserve the accuracy of time recovery
        
        % BPF
        tdelay = 1/(4*sim.Rs); % 90deg delay at symbol rate
        clk = filtfilt(TimeRec.bpf.num, TimeRec.bpf.den, ynl); % use zero-phase filtering
        clk = real(ifft(fft(clk).*ifftshift(exp(-1j*2*pi*f*tdelay)))); % 90deg shift at symbol rate   
        
        % Limiting amps
        clk(clk > 0) = 1;
        clk(clk < 0) = 0;
        
        % Rising-edge detection
        samp = diff([clk 0]);
        idx = find(samp > 0);
        
        % discard adaptation period
        idx(diff([idx 0]) < round(sim.Mct*0.4)) = [];
                
        % Sampling
        tsamp = t(idx);
        ysamp = trim_pad(yct(idx));
        
        Ndiscard = TimeRec.Ndiscard;
        
        % Plots
        if exist('verbose', 'var') && verbose
            try % sometimes matrix length will not match
                figure(405), clf, hold on
                plot(1e12*(sim.t(1:sim.Mct:min(length(sim.t), length(tsamp)*sim.Mct)) - tsamp))
                xlabel('Symbol')
                ylabel('Sampling time error (ps)')
                title('Timing recovery')
                drawnow
            catch e
                disp('Error while trying to plot Sampling time error')
            end
        end
 
    case 'spectral-line-pll'
    %% Spectral line syncronizer based on squaring and PLL
        % Differentiator -> nonlinearity (squaring) -> band pass filter
        % centered at the symbol rate
        % Note: BPF could be replaced by a PLL, but this is not considered
        % here        
        
        % Upsample signal to have better time accuray
        yct = resample(yct, TimeRec.Mct, sim.Mct);
        fs = sim.fs*TimeRec.Mct/sim.Mct;
        [~, tct] = freq_time(length(yct), fs);
        
        Squarer = AnalogSquaring(TimeRec.squarerFilt, TimeRec.N0, fs); % Squarer
        M = AnalogMixer(TimeRec.mixerFilt, TimeRec.N0, fs); % Mixer
        
        % Loop filter
        nums = [2*TimeRec.csi*TimeRec.wn TimeRec.wn^2];
        dens = [1 0 0]; % descending powers of s
        if strcmpi(TimeRec.CT2DT, 'bilinear')
            [numz, denz] = bilinear(nums, dens, fs); 
        elseif strcmpi(TimeRec.CT2DT, 'impinvar')
            [numz, denz] = impinvar(nums, dens, fs); 
        else
            error('analog_time_recovery: continuous-time to discrete-time conversion method (CT2DT) must be either bilinear or impinvar')
        end         
        LoopFilter = ClassFilter(numz, denz, fs);
        
        % Only uses I component in estimating clock
        y = real(yct);
              
        % Nonlinearity
        ynl = Squarer.square(y);

        % Note: group delay due to filtering operations is not removed
        f0 = sim.Rs; % VCO natural frequency
        S = zeros(size(tct));
        Sf = zeros(size(tct));
        clk = zeros(size(tct));
        additionalDelay = max(1, TimeRec.loopDelay); % at least 1 sample
        fprintf('analog_group_delay: Time recovery loop delay = %.2f ps (%d samples)\n', additionalDelay/fs*1e12, additionalDelay);
        for t = additionalDelay+1:length(tct)
            % VCO: generates VCO output
            clk(t) = sin(2*pi*f0*tct(t) - Sf(t-additionalDelay));

            % Mixing            
            S(t) = M.mix(ynl(t), clk(t));
            
            % Loop filter
            Sf(t) = LoopFilter.filter(S(t));  
        end
              
        % Limiting amps
        clk(clk > 0) = 1;
        clk(clk < 0) = 0;
        
        % Rising-edge detection
        samp = diff([clk 0]);
        idx = find(samp > 0);
        
        % Sampling
        tsamp = tct(idx);
        ysamp = trim_pad(yct(idx));   
        
        Ndiscard = TimeRec.Ndiscard;
        
        % Plots
        if exist('verbose', 'var') && verbose
            try % sometimes matrix length will not match
                figure(405), clf
                subplot(211), hold on, box on
                tplot = sim.t(1:sim.Mct:min(length(sim.t), length(tsamp)*sim.Mct));
                plot(1e12*(tplot - tsamp(1:min(length(tsamp), length(tplot)))))
                a = axis;
                plot((sim.Ndiscard+Ndiscard(1)+1)*[1 1], a(3:4), ':k')
                plot((sim.Nsymb-sim.Ndiscard-Ndiscard(2))*[1 1], a(3:4), ':k')
                xlabel('Symbol')
                ylabel('Sampling time error (ps)')
                title('Timing recovery')
                axis tight
                subplot(212)
                plot(tct*1e9, Sf)
                xlabel('Time (ns)')
                ylabel('Loop filter output')
                drawnow
            catch e
                disp('Error while trying to plot Sampling time error')
            end
        end
 
    case 'none'
        idx = 1:sim.Mct:length(yct);
        ysamp = yct(idx);
        Ndiscard = [0 0];
    otherwise
        error('analog_time_recovery: invalide time recovery type')
end

function y = trim_pad(x)
%% Ensures that vector x has sim.Nsymb points by trimming or padding zeros
    N = sim.Nsymb;
    if length(x) > N
        fprintf('analog_time_recovery: first %d samples were discarded\n', length(x)-N)
        y = x(length(x)-N+1:end);
    elseif length(x) < N
        fprintf('analog_time_recovery: %d zeros were padded at the end\n', N-length(x))
        y = [x zeros(1, N-length(x))];
    else
        y = x;
    end
end
end
