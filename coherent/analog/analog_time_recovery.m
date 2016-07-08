function [ysamp, idx] = analog_time_recovery(yct, TimeRec, sim)
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

% Only I component is used in timing recovery
y = real(yct);

switch lower(TimeRec.type)
    case 'squarer'
        %% Spectral line syncronizer based on squaring and bandpass filtering
        % Differentiator -> nonlinearity (squaring) -> band pass filter
        % centered at the symbol rate
        % Note: BPF could be replaced by a PLL, but this is not considered
        % here
        
        ynl = diff([y 0]).^2;
        clk = real(ifft(fft(ynl).*ifftshift(TimeRec.bpf.H(sim.f/sim.fs))));
        clk(clk > 0) = 1;
        clk(clk < 0) = 0;
        samp = diff([clk 0]);
        idx = find(samp > 0);
        ysamp = yct(idx);
                
    case 'mueller-muller'
        %% Mueller and Muller (M&M) algorithm
        % Discrete-time error-tracking algorithm     
        % Build loop filter
        csi = TimeRec.csi;
        wn = TimeRec.wn;
        if not(isfield(TimeRec, 'CT2DT'))
            TimeRec.CT2DT = 'bilinear';
        end
        
        % Open-loop analog filter coefficients
        nums = [0 2*csi*wn wn^2];
        dens = [1 0 0]; % descending powers of s
        % Closed-loop digital filter coefficients using bilinear transformation
        if strcmpi(TimeRec.CT2DT, 'bilinear')
            [numz, denz] = bilinear(nums, dens+nums, sim.Rs); % ascending powers of z^–1
        elseif strcmpi(TimeRec.CT2DT, 'impinvar')
            [numz, denz] = impinvar(nums, dens+nums, sim.Rs); % ascending powers of z^–1
        else
            error('analog_time_recovery: invalid method for continuous time to discrete time conversion')
        end
        
        LoopFilter = ClassFilter(numz, denz, sim.Rs);
        
        T = sim.Mct;
        xhat = [0 0];
        k = 2;
        n = T+1;
        s = zeros(1, sim.Nsymb);
        idx = ones(1, sim.Nsymb);
        mk = zeros(1, sim.Nsymb);
        ysamp = zeros(1, sim.Nsymb);
        while n <= length(y) && k <= sim.Nsymb
            ysamp(k) = y(n);
            idx(k) = n;
            xhat(2) = xhat(1);
            xhat(1) = TimeRec.detect(y(n)); % make symbol decision
            mk(k) = xhat(2)*ysamp(k) - xhat(1)*ysamp(k-1); % M&M timing error detector

            % Loop filter 
            s(k) = LoopFilter.filter(mk(k));

            s(k) = max(-T+1, s(k)); % makes sure sampler is always moving forward
            n = round(n + T + s(k));
            k = k + 1;
        end

        ysamp = yct(idx(idx ~= 0));
       
    case 'none'
        if mod(sim.Mct, 2) == 0 % even
            idx = sim.Mct/2:sim.Mct:length(y);
        else % odd
            idx = (sim.Mct + 1)/2:sim.Mct:length(y);
        end
        
        ysamp = yct(idx);     
    otherwise
        error('analog_time_recovery: invalide time recovery type')
end